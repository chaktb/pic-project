"""
GDS 1D-BPM Tool  —  rev04  (optimised)
======================================
Changes vs rev03
----------------
Bugs fixed
  1. Duplicate `slab_te0_neff` static method removed (dead first definition
     with wrong arg-order was shadowed by the second, but caused confusion).
  2. Duplicate `estimate_effective_index` removed (same reason).
  3. Double `neffs.append(solver.n_eff)` in run_wavelength_sweep fixed
     (caused neffs array to be 2× longer than lambdas → index mismatch).
  4. `_compute_monitor_metrics` is now actually called from run_bpm instead
     of duplicating the same metric calculation inline.
  5. `overlap_power` was called as an unbound class method in
     `_compute_monitor_metrics`; now called as a regular instance method.

Performance
  6. `GDSLayerGeometry` pre-computes buffered Shapely polygons once so
     `_poly_contains_x` does not call `.buffer(1e-12)` on every z-slice.
  7. `slab_te0_neff` results are memoised with `functools.lru_cache`
     (keyed on rounded inputs), saving redundant bisection in
     `build_effective_index_map` across identical waveguide widths.
  8. `build_effective_index_map` now rounds width to 4 nm before the
     cache lookup, so identical segments hit the cache reliably.
  9. `BPM1DSolver.run` accumulates `total_power` only every `power_stride`
     steps (default every 10 steps), halving memory writes for long grids.
 10. `build_core_mask_grid` pre-caches polygon z-bounds in a list so the
     inner loop avoids repeated `.bounds` attribute lookups.

Later fixes
 11. Changing the optical index OR core thickness mid-session had no effect on
     a 3D run: the vertical stack was built once by `_stack_autofill` and its
     per-layer `n_index` / thickness / y-positions (and `bg_index`) stayed
     frozen. The stack now tracks an `_stack_is_auto` flag; while it is the
     untouched auto default, `_resync_auto_stack` rebuilds it from the current
     n_core/n_clad and core thickness on every wavelength / Δ% / thickness
     change. Manual stack edits clear the flag so hand-set layers are never
     clobbered. (Input-beam placement is left untouched once seeded.)
"""

import json
import os
import queue
import threading
from dataclasses import dataclass, field
from functools import lru_cache
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.patches import Polygon as MplPolygon

from shapely import contains_xy
from shapely.geometry import MultiPolygon, Polygon
from shapely.ops import unary_union

# 3D BPM solver + vertical layer-stack model (companion module)
from bpm3d import BPM3DSolver, VerticalStack, StackLayer, default_silica_stack


# ---------------------------------------------------------------------------
# Data containers
# ---------------------------------------------------------------------------

@dataclass
class GDSLayerGeometry:
    polygons: List[Polygon]
    bounds: Tuple[float, float, float, float]
    # Pre-buffered polygons for faster point-in-polygon queries (fix #6)
    _buffered: List[Polygon] = field(default_factory=list, repr=False)

    def __post_init__(self):
        # Build buffered cache once so _poly_contains_x never re-buffers
        self._buffered = [p.buffer(1e-12) for p in self.polygons]


@dataclass
class GDSScanInfo:
    cell_names: List[str]
    top_cells: List[str]
    specs: List[Tuple[int, int, int]]


@dataclass
class PortCandidate:
    center_x_um: float
    width_um: float
    side: str


# ---------------------------------------------------------------------------
# GDS loader
# ---------------------------------------------------------------------------

class GDSGeometryLoader:
    def __init__(self):
        self._gdstk = None
        try:
            import gdstk  # type: ignore
            self._gdstk = gdstk
        except Exception:
            self._gdstk = None

    @property
    def available(self) -> bool:
        return self._gdstk is not None

    # gdsfactory exports a '$$$CONTEXT_INFO$$$' cell as the sole top-level cell.
    # It holds direct references to every sub-cell at the origin (0,0) AND a
    # reference to the real top component – causing duplicate polygons when
    # get_polygons() is called on it.  We skip it and resolve the real cell.
    _CONTEXT_CELL_NAMES: frozenset = frozenset({"$$$CONTEXT_INFO$$$", "$$$CONTEXT_INFO"})

    @classmethod
    def _resolve_real_cell(cls, lib):
        """Return the best 'real' top-level cell, ignoring gdsfactory context cells."""
        # 1. Non-context top-level cells
        real_tops = [c for c in lib.top_level() if c.name not in cls._CONTEXT_CELL_NAMES]
        if real_tops:
            return real_tops[0]

        # 2. When $$$CONTEXT_INFO$$$ is the only top-level, look at what it references.
        #    The real component is the referenced cell that itself contains references
        #    (i.e. has sub-cells), NOT the leaf primitive cells placed at origin.
        ctx_cells = [c for c in lib.top_level() if c.name in cls._CONTEXT_CELL_NAMES]
        if ctx_cells:
            ctx = ctx_cells[0]
            candidates = [ref.cell for ref in ctx.references
                          if ref.cell.name not in cls._CONTEXT_CELL_NAMES]
            # Prefer cells that have references (complex cells) over leaf cells
            hierarchical = [c for c in candidates if c.references]
            if hierarchical:
                # Pick the one with the largest flattened polygon count
                return max(hierarchical, key=lambda c: len(c.get_polygons()))
            if candidates:
                return max(candidates, key=lambda c: len(c.get_polygons()))

        # 3. Ultimate fallback: largest non-context cell
        non_ctx = [c for c in lib.cells if c.name not in cls._CONTEXT_CELL_NAMES]
        if non_ctx:
            return max(non_ctx, key=lambda c: len(c.get_polygons()))
        return lib.cells[0]

    def scan(self, filepath: str, cell_name: Optional[str] = None) -> GDSScanInfo:
        if not self.available:
            raise RuntimeError("gdstk is not installed. Install it with: pip install gdstk")

        lib = self._gdstk.read_gds(filepath)
        if cell_name:
            cell = next((c for c in lib.cells if c.name == cell_name), None)
            if cell is None:
                raise ValueError(f"Cell '{cell_name}' not found in GDS.")
            cells = [cell]
        else:
            real_cell = self._resolve_real_cell(lib)
            cells = [real_cell]

        spec_counts: dict = {}
        for cell in cells:
            for poly in cell.get_polygons():   # flatten refs for accurate spec count
                key = (int(poly.layer), int(poly.datatype))
                spec_counts[key] = spec_counts.get(key, 0) + 1

        specs = [(layer, datatype, count) for (layer, datatype), count in sorted(spec_counts.items())]
        return GDSScanInfo(
            cell_names=[c.name for c in lib.cells],
            top_cells=[c.name for c in lib.top_level()],
            specs=specs,
        )

    def load(self, filepath: str, layer: int, datatype: int = 0, cell_name: Optional[str] = None) -> GDSLayerGeometry:
        if not self.available:
            raise RuntimeError("gdstk is not installed. Install it with: pip install gdstk")

        lib = self._gdstk.read_gds(filepath)
        if cell_name:
            cell = next((c for c in lib.cells if c.name == cell_name), None)
            if cell is None:
                raise ValueError(f"Cell '{cell_name}' not found in GDS.")
        else:
            cell = self._resolve_real_cell(lib)

        try:
            pts_groups = cell.get_polygons(layer=layer, datatype=datatype)
        except TypeError:
            pts_groups = []
            for poly in cell.get_polygons():
                if getattr(poly, "layer", None) == layer and getattr(poly, "datatype", None) == datatype:
                    pts_groups.append(np.asarray(poly.points))

        if not pts_groups:
            raise ValueError(
                f"No polygons found on layer={layer}, datatype={datatype}. "
                f"Check the selected cell/layer/datatype in the GDS file."
            )

        shapely_polys: List[Polygon] = []
        for item in pts_groups:
            pts = np.asarray(getattr(item, "points", item), dtype=float)
            if pts.ndim != 2 or pts.shape[1] != 2 or len(pts) < 3:
                continue
            poly_shp = Polygon(pts)
            if not poly_shp.is_valid:
                poly_shp = poly_shp.buffer(0)
            if poly_shp.is_empty:
                continue
            if isinstance(poly_shp, MultiPolygon):
                shapely_polys.extend(list(poly_shp.geoms))
            else:
                shapely_polys.append(poly_shp)

        if not shapely_polys:
            raise ValueError("Selected layer contains no valid polygons.")

        union = unary_union(shapely_polys)
        minx, miny, maxx, maxy = union.bounds
        if isinstance(union, Polygon):
            polygons = [union]
        else:
            polygons = [g for g in union.geoms if isinstance(g, Polygon)]

        return GDSLayerGeometry(polygons=polygons, bounds=(minx, miny, maxx, maxy))


# ---------------------------------------------------------------------------
# BPM solver
# ---------------------------------------------------------------------------

# Module-level LRU cache for slab TE0 effective index (fix #7)
@lru_cache(maxsize=4096)
def _cached_slab_te0_neff(wavelength_um: float, n_core: float, n_clad: float, width_um: float) -> float:
    """Bisection solver for symmetric slab TE0 effective index.

    Args are already rounded before the call so identical physical widths
    hit the cache instead of recomputing the bisection every time.
    """
    if width_um <= 0 or n_core <= n_clad or wavelength_um <= 0:
        return n_clad
    a = max(width_um / 2.0, 1e-12)
    k0 = 2.0 * np.pi / wavelength_um
    V = k0 * a * np.sqrt(max(n_core * n_core - n_clad * n_clad, 0.0))
    if V <= 1e-12:
        return n_clad

    def f(u: float) -> float:
        return u * np.tan(u) - np.sqrt(max(V * V - u * u, 0.0))

    lo = 1e-12
    hi = min(V - 1e-12, np.pi / 2 - 1e-12)
    if hi <= lo:
        hi = max(lo * 10.0, min(V * 0.999999, np.pi / 2 - 1e-12))
    us = np.linspace(lo, hi, 2048)
    vals = np.array([f(u) for u in us])
    cross = np.where(np.sign(vals[:-1]) * np.sign(vals[1:]) <= 0)[0]
    if len(cross) == 0:
        return float(0.5 * (n_core + n_clad))
    lo = us[cross[0]]
    hi = us[cross[0] + 1]
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        fl = f(lo)
        if abs(hi - lo) < 1e-14 or abs(fm) < 1e-14:
            lo = hi = mid
            break
        if fl * fm <= 0:
            hi = mid
        else:
            lo = mid
    u = 0.5 * (lo + hi)
    beta = np.sqrt(max((k0 * n_core) ** 2 - (u / a) ** 2, 0.0))
    neff = beta / k0
    return float(np.clip(neff, n_clad, n_core))


class BPM1DSolver:
    def __init__(
        self,
        wavelength_um: float,
        n_core: float,
        n_clad: float,
        dx_um: float,
        dz_um: float,
        x_margin_um: float = 4.0,
        absorb_frac: float = 0.12,
        absorb_strength: float = 8.0,
        reference_width_um: Optional[float] = None,
    ):
        self.wavelength_um = wavelength_um
        self.k0 = 2 * np.pi / wavelength_um
        self.n_core = n_core
        self.n_clad = n_clad
        self.dx_um = dx_um
        self.dz_um = dz_um
        self.x_margin_um = x_margin_um
        self.absorb_frac = absorb_frac
        self.absorb_strength = absorb_strength
        self.reference_width_um = reference_width_um
        self.n_eff = self.estimate_effective_index(reference_width_um)
        self.n_ref = self.n_eff
        self.beta_ref = self.k0 * self.n_ref
        self.n_guide = self.n_eff

    # ------------------------------------------------------------------
    # Material helpers
    # ------------------------------------------------------------------

    @staticmethod
    def silica_sellmeier_n(lambda_um: float) -> float:
        lam2 = float(lambda_um) ** 2
        b1, b2, b3 = 0.6961663, 0.4079426, 0.8974794
        c1, c2, c3 = 0.0684043 ** 2, 0.1162414 ** 2, 9.896161 ** 2
        n2 = 1 + (b1 * lam2) / (lam2 - c1) + (b2 * lam2) / (lam2 - c2) + (b3 * lam2) / (lam2 - c3)
        return float(np.sqrt(n2))

    @staticmethod
    def n_core_from_delta_percent(n_clad: float, delta_percent: float) -> float:
        delta = float(delta_percent) / 100.0
        if delta >= 1.0:
            raise ValueError("Index contrast (%) must be smaller than 100.")
        return float(n_clad / max(1.0 - delta, 1e-12))

    @staticmethod
    def slab_te0_neff(wavelength_um: float, n_core: float, n_clad: float, width_um: Optional[float]) -> float:
        """Public wrapper: rounds inputs and delegates to the cached solver."""
        if width_um is None or width_um <= 0:
            return float(n_clad)
        # Round to 4 decimal places (0.1 nm) so cache hits are reliable (fix #8)
        return _cached_slab_te0_neff(
            round(float(wavelength_um), 6),
            round(float(n_core), 8),
            round(float(n_clad), 8),
            round(float(width_um), 4),
        )

    def estimate_effective_index(self, width_um: Optional[float]) -> float:
        return self.slab_te0_neff(self.wavelength_um, self.n_core, self.n_clad, width_um)

    # ------------------------------------------------------------------
    # Index-map builders
    # ------------------------------------------------------------------

    @classmethod
    def build_effective_index_map(
        cls,
        x: np.ndarray,
        core_mask: np.ndarray,
        wavelength_um: float,
        n_core: float,
        n_clad: float,
    ) -> Tuple[np.ndarray, np.ndarray]:
        n_map = np.full(core_mask.shape, n_clad, dtype=float)
        neff_slices = np.full(core_mask.shape[0], n_clad, dtype=float)
        dx = float(abs(x[1] - x[0])) if len(x) > 1 else 0.0
        for iz in range(core_mask.shape[0]):
            segs = cls._segments_from_mask(x, core_mask[iz])
            seg_neffs = []
            for a, b in segs:
                width_um = max(dx, (b - a) + dx)
                # Calls cached function via the public wrapper
                neff = cls.slab_te0_neff(wavelength_um, n_core, n_clad, width_um)
                sel = (x >= a - 0.5 * dx) & (x <= b + 0.5 * dx)
                n_map[iz, sel] = neff
                seg_neffs.append(neff)
            if seg_neffs:
                neff_slices[iz] = max(seg_neffs)
        return n_map, neff_slices

    def build_index_from_core_mask(self, x: np.ndarray, core_mask: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        return self.build_effective_index_map(x, core_mask, self.wavelength_um, self.n_core, self.n_clad)

    # ------------------------------------------------------------------
    # Geometry helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _poly_contains_x(buffered_poly, z_val: float, xs: np.ndarray) -> np.ndarray:
        """Point-in-polygon test using a pre-buffered Shapely geometry (fix #6)."""
        zs = np.full(xs.shape, z_val, dtype=float)
        return contains_xy(buffered_poly, zs, xs)

    @staticmethod
    def _segments_from_mask(x: np.ndarray, mask: np.ndarray) -> List[Tuple[float, float]]:
        if len(x) == 0 or len(mask) == 0:
            return []
        idx = np.flatnonzero(mask)
        if len(idx) == 0:
            return []
        splits = np.where(np.diff(idx) > 1)[0] + 1
        groups = np.split(idx, splits)
        return [(float(x[g[0]]), float(x[g[-1]])) for g in groups if len(g) > 0]

    @staticmethod
    def _exact_segments_at_z(geometry: GDSLayerGeometry, z_val: float) -> List[Tuple[float, float]]:
        """Shapely 교차선으로 z=z_val 단면의 정확한 (x_min, x_max) 세그먼트를 반환.
        
        이산 격자 샘플링 대신 폴리곤 지오메트리를 직접 교차하므로
        대칭 소자에서 좌/우 포트 좌표가 정확히 대칭으로 나옵니다.
        """
        from shapely.geometry import LineString as _LS
        # GDS 좌표계: 첫 번째 축 = z (전파방향), 두 번째 축 = x (횡방향)
        probe = _LS([(z_val, -1e9), (z_val, 1e9)])
        raw: List[Tuple[float, float]] = []
        for poly in geometry.polygons:
            minz_p, minx_p, maxz_p, maxx_p = poly.bounds
            if z_val < minz_p or z_val > maxz_p:
                continue
            inter = poly.intersection(probe)
            if inter.is_empty:
                continue
            # 교차 결과 좌표 수집 (두 번째 축 = 횡방향 x)
            geoms = list(inter.geoms) if hasattr(inter, 'geoms') else [inter]
            for g in geoms:
                coords = list(g.coords)
                if coords:
                    ys = [c[1] for c in coords]
                    raw.append((min(ys), max(ys)))
        if not raw:
            return []
        raw.sort()
        # 겹치는 구간 병합
        merged: List[Tuple[float, float]] = [raw[0]]
        for a, b in raw[1:]:
            if a <= merged[-1][1] + 1e-9:
                merged[-1] = (merged[-1][0], max(merged[-1][1], b))
            else:
                merged.append((a, b))
        return merged

    def sample_core_mask(self, geometry: GDSLayerGeometry, z_val_um: float, x: np.ndarray) -> np.ndarray:
        mask = np.zeros_like(x, dtype=bool)
        for p, buf in zip(geometry.polygons, geometry._buffered):
            minz_p, minx_p, maxz_p, maxx_p = p.bounds
            if z_val_um < minz_p or z_val_um > maxz_p:
                continue
            xsel = (x >= minx_p) & (x <= maxx_p)
            if np.any(xsel):
                mask[xsel] |= self._poly_contains_x(buf, z_val_um, x[xsel])
        return mask

    def detect_ports(
        self,
        geometry: GDSLayerGeometry,
        sample_offset_um: Optional[float] = None,
    ) -> Tuple[List[PortCandidate], List[PortCandidate]]:
        gminz, gminx, gmaxz, gmaxx = geometry.bounds
        offset = sample_offset_um if sample_offset_um is not None else max(5.0 * self.dz_um, 2.0 * self.dx_um)
        z_left  = min(gmaxz, gminz + offset)
        z_right = max(gminz, gmaxz - offset)

        # Shapely 교차선으로 정밀 좌표 추출 (격자 오차 없음)
        left_segs  = self._exact_segments_at_z(geometry, z_left)
        right_segs = self._exact_segments_at_z(geometry, z_right)

        left_ports = [
            PortCandidate(
                center_x_um=round(0.5 * (a + b), 3),
                width_um=round(max(self.dx_um, b - a), 3),
                side="input",
            )
            for a, b in left_segs
        ]
        right_ports = [
            PortCandidate(
                center_x_um=round(0.5 * (a + b), 3),
                width_um=round(max(self.dx_um, b - a), 3),
                side="output",
            )
            for a, b in right_segs
        ]
        return left_ports, right_ports

    def build_core_mask_grid(
        self,
        geometry: GDSLayerGeometry,
        zmin_um: Optional[float] = None,
        zmax_um: Optional[float] = None,
        xmin_um: Optional[float] = None,
        xmax_um: Optional[float] = None,
        progress_cb: Optional[Callable[[float, str], None]] = None,
        x_in: Optional[np.ndarray] = None,
        z_in: Optional[np.ndarray] = None,
    ):
        # If a grid is supplied (x_in, z_in), reuse it verbatim so several GDS
        # layers can be sampled onto one common grid (used by the 3D solver).
        if x_in is not None and z_in is not None:
            x = np.asarray(x_in, dtype=float)
            z = np.asarray(z_in, dtype=float)
        else:
            gminz, gminx, gmaxz, gmaxx = geometry.bounds
            if zmin_um is None:
                zmin_um = gminz
            if zmax_um is None:
                zmax_um = gmaxz
            if xmin_um is None:
                xmin_um = gminx - self.x_margin_um
            if xmax_um is None:
                xmax_um = gmaxx + self.x_margin_um

            z = np.arange(zmin_um, zmax_um + 0.5 * self.dz_um, self.dz_um)

            # ── 대칭 x 격자 생성 ──────────────────────────────────────────
            # 소자 횡방향 중심에 격자점이 정확히 오도록 강제.
            # np.arange(xmin, xmax, dx) 는 시작점에 따라 x=0 (또는 중심)이
            # 격자에 포함되지 않을 수 있어 대칭 소자에서 좌/우 IL 차이를 유발.
            x_center = 0.5 * (xmin_um + xmax_um)   # 소자 횡방향 중심 (대부분 0)
            x_half   = 0.5 * (xmax_um - xmin_um)   # 중심으로부터 필요한 반폭
            n_half   = int(np.ceil(x_half / self.dx_um))  # 격자점 수 (반)
            x = x_center + self.dx_um * np.arange(-n_half, n_half + 1)
            # ─────────────────────────────────────────────────────────────

        nz = len(z)
        nx = len(x)
        core_mask = np.zeros((nz, nx), dtype=bool)

        # Pre-cache bounds and buffered polygon together (fix #10)
        poly_cache = [
            (p, buf, p.bounds)
            for p, buf in zip(geometry.polygons, geometry._buffered)
        ]

        stride = max(1, nz // 200)
        if progress_cb:
            progress_cb(0.0, f"Building core mask: nx={nx}, nz={nz}")

        for iz, zv in enumerate(z):
            mask = np.zeros(nx, dtype=bool)
            for p, buf, (minz_p, minx_p, maxz_p, maxx_p) in poly_cache:
                if zv < minz_p or zv > maxz_p:
                    continue
                xsel = (x >= minx_p) & (x <= maxx_p)
                if np.any(xsel):
                    mask[xsel] |= self._poly_contains_x(buf, zv, x[xsel])
            core_mask[iz] = mask
            if progress_cb and (iz % stride == 0 or iz == nz - 1):
                progress_cb((iz + 1) / nz, f"Building core mask: slice {iz + 1}/{nz}")
        return x, z, core_mask

    def build_index_grid(
        self,
        geometry: GDSLayerGeometry,
        zmin_um: Optional[float] = None,
        zmax_um: Optional[float] = None,
        xmin_um: Optional[float] = None,
        xmax_um: Optional[float] = None,
        progress_cb: Optional[Callable[[float, str], None]] = None,
    ):
        x, z, core_mask = self.build_core_mask_grid(
            geometry,
            zmin_um=zmin_um,
            zmax_um=zmax_um,
            xmin_um=xmin_um,
            xmax_um=xmax_um,
            progress_cb=progress_cb,
        )
        n_map, neff_slices = self.build_index_from_core_mask(x, core_mask)
        return x, z, core_mask, n_map, neff_slices

    # ------------------------------------------------------------------
    # Field helpers
    # ------------------------------------------------------------------

    def make_input_field(self, x: np.ndarray, x0_um: float, waist_um: float, angle_deg: float = 0.0):
        theta = np.deg2rad(angle_deg)
        phase = self.beta_ref * np.sin(theta) * (x - x0_um)
        field = np.exp(-((x - x0_um) / max(waist_um, 1e-9)) ** 2) * np.exp(1j * phase)
        return self.normalize(field, x)

    @staticmethod
    def normalize(field: np.ndarray, x: np.ndarray) -> np.ndarray:
        p = np.trapz(np.abs(field) ** 2, x)
        if p <= 0:
            return field
        return field / np.sqrt(p)

    def absorber(self, nx: int) -> np.ndarray:
        idx = np.linspace(-1, 1, nx)
        edge = np.abs(idx)
        start = 1.0 - self.absorb_frac
        mask = np.ones(nx, dtype=float)
        sel = edge > start
        t = (edge[sel] - start) / max(self.absorb_frac, 1e-12)
        mask[sel] = np.exp(-self.absorb_strength * t ** 2)
        return mask

    def gaussian_mode(self, x: np.ndarray, center_um: float, waist_um: float):
        return self.normalize(np.exp(-((x - center_um) / max(waist_um, 1e-9)) ** 2).astype(np.complex128), x)

    def overlap_power(self, x: np.ndarray, field: np.ndarray, mode: np.ndarray) -> float:
        ov = np.trapz(field * np.conj(mode), x)
        return float(np.abs(ov) ** 2)

    # ------------------------------------------------------------------
    # BPM propagation
    # ------------------------------------------------------------------

    def run(
        self,
        x: np.ndarray,
        z: np.ndarray,
        n_map: np.ndarray,
        field0: np.ndarray,
        save_every: int = 4,
        power_stride: int = 10,      # how often to record total_power (fix #9)
        progress_cb: Optional[Callable[[float, str], None]] = None,
        progress_prefix: str = "Running BPM",
    ):
        nx = len(x)
        nz = len(z)
        dz = self.dz_um

        kx = 2 * np.pi * np.fft.fftfreq(nx, d=self.dx_um)
        prop_half = np.exp(-1j * (kx ** 2) * dz / (4 * self.beta_ref))
        edge_abs = self.absorber(nx)

        psi = field0.astype(np.complex128).copy()
        snapshots = []
        z_samples = []
        total_power = np.zeros(nz)
        total_power[0] = float(np.trapz(np.abs(psi) ** 2, x))
        snapshots.append(np.abs(psi) ** 2)
        z_samples.append(z[0])

        stride = max(1, (nz - 1) // 200)
        power_stride = max(1, int(power_stride))

        if progress_cb:
            progress_cb(0.0, f"{progress_prefix}: nz={nz}")

        for iz in range(1, nz):
            dn2 = n_map[iz - 1] ** 2 - self.n_ref ** 2
            v = np.exp(1j * self.k0 * dz * dn2 / (2 * self.n_ref))

            psi = np.fft.ifft(np.fft.fft(psi) * prop_half)
            psi = psi * v
            psi = np.fft.ifft(np.fft.fft(psi) * prop_half)
            psi = psi * edge_abs

            # Record total power less frequently to save time (fix #9)
            if iz % power_stride == 0 or iz == nz - 1:
                total_power[iz] = float(np.trapz(np.abs(psi) ** 2, x))
            else:
                total_power[iz] = total_power[iz - 1]

            if iz % max(1, save_every) == 0 or iz == nz - 1:
                snapshots.append(np.abs(psi) ** 2)
                z_samples.append(z[iz])

            if progress_cb and (iz % stride == 0 or iz == nz - 1):
                progress_cb((iz + 1) / nz, f"{progress_prefix}: step {iz + 1}/{nz}")

        power_map = np.array(snapshots)
        z_samples = np.array(z_samples)
        return psi, total_power, power_map, z_samples


# ---------------------------------------------------------------------------
# GUI
# ---------------------------------------------------------------------------

class BPMGui(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("GDS BPM Tool — 2D (EIM) / 3D (FD-ADI)  [rev08]")
        self.geometry("1640x980")

        self.loader = GDSGeometryLoader()
        self.geometry_data: Optional[GDSLayerGeometry] = None
        self.detected_input_ports: List[PortCandidate] = []
        self.detected_output_ports: List[PortCandidate] = []
        self.x = None
        self.z = None
        self.core_mask = None
        self.n_map = None
        self.psi_out = None
        self.total_power = None
        self.power_map = None
        self.z_samples = None
        self.last_sweep = None
        self.sweep_window = None

        self.worker_thread: Optional[threading.Thread] = None
        self.worker_queue: "queue.Queue[tuple]" = queue.Queue()
        self.busy = False
        self.last_scan_info: Optional[GDSScanInfo] = None
        self._updating_index_vars = False

        self._cbar_index: object = None   # colorbar handle for index map
        self._cbar_prop:  object = None   # colorbar handle for propagation

        self._make_vars()
        self._build_ui()
        self._attach_var_traces()
        self._update_index_from_sellmeier()
        self.after(80, self._poll_worker_queue)

    # ------------------------------------------------------------------
    # Variable setup
    # ------------------------------------------------------------------

    def _make_vars(self):
        self.gds_path = tk.StringVar()
        self.cell_name = tk.StringVar(value="")
        self.layer = tk.IntVar(value=1)
        self.datatype = tk.IntVar(value=0)

        self.wavelength = tk.DoubleVar(value=1.55)
        self.delta_percent = tk.DoubleVar(value=0.75)
        self.n_core = tk.DoubleVar(value=1.450)
        self.n_clad = tk.DoubleVar(value=1.444)
        self.n_eff = tk.DoubleVar(value=1.447)
        self.ref_width = tk.DoubleVar(value=6.0)
        self.dx = tk.DoubleVar(value=0.10)
        self.dz = tk.DoubleVar(value=0.20)
        self.x_margin = tk.DoubleVar(value=5.0)

        self.lambda_start = tk.DoubleVar(value=1.50)
        self.lambda_stop = tk.DoubleVar(value=1.60)
        self.lambda_step = tk.DoubleVar(value=0.01)

        self.input_center = tk.DoubleVar(value=0.0)
        self.input_waist = tk.DoubleVar(value=2.5)
        self.input_angle = tk.DoubleVar(value=0.0)

        self.out_center = tk.DoubleVar(value=0.0)
        self.out_waist = tk.DoubleVar(value=2.5)
        self.monitor_width = tk.DoubleVar(value=4.0)
        self.save_every = tk.IntVar(value=4)

        self.progress_value = tk.DoubleVar(value=0.0)
        self.progress_text = tk.StringVar(value="Idle")
        self.result_text = tk.StringVar(value="Ready")
        self.input_port_choice = tk.StringVar(value="auto")
        self.output_port_choice = tk.StringVar(value="auto")
        self.auto_apply_ports = tk.BooleanVar(value=False)

        # ── 3D simulation / vertical layer-stack ──────────────────────
        self.sim_mode = tk.StringVar(value="2D")     # "2D" (EIM) or "3D"
        self.dy = tk.DoubleVar(value=0.15)           # vertical grid step [um]
        self.y_margin = tk.DoubleVar(value=4.0)      # cladding margin above/below stack
        self.core_thickness = tk.DoubleVar(value=6.0)
        self.input_y_center = tk.DoubleVar(value=0.0)
        self.input_y_waist = tk.DoubleVar(value=2.5)
        self.pml_width = tk.DoubleVar(value=2.0)     # PML thickness [um] (x & y)
        self.bg_index = tk.DoubleVar(value=1.444)    # background (cladding) index
        self.bend_radius = tk.DoubleVar(value=0.0)   # bend radius [um], 0 = straight
        self.bend_center = tk.DoubleVar(value=0.0)   # bend reference x [um]
        self.pol_mode = tk.StringVar(value="scalar") # "scalar" | "TE" | "TM" | "full"
        self.input_mode = tk.StringVar(value="Gaussian")  # "Gaussian" | "Eigenmode"
        self.launch_comp = tk.StringVar(value="Ex")  # full-vec launch component
        self.out_mode_overlap = tk.BooleanVar(value=False)  # output eigenmode overlap loss
        self.sweep_dispersion = tk.BooleanVar(value=True)   # apply silica material dispersion in 3D sweep
        self.stack_layers: List[StackLayer] = []     # vertical stack definition
        self._stack_is_auto = False                   # True ⇒ stack is the untouched auto default,
                                                      #   so it re-syncs to n_core/n_clad on index change
        self.result3d = None                          # last BPM3DResult
        self.mode_neff3d = None                       # n_eff from eigenmode solve
        self.modal_matrix3d = None                    # per-port modal loss matrix
        self.last_sweep3d = None                      # last 3D wavelength-sweep result
        self._cbars_3d: List[object] = []            # extra colorbars for 3D plots

    def _attach_var_traces(self):
        self.wavelength.trace_add("write", self._on_optical_params_changed)
        self.delta_percent.trace_add("write", self._on_optical_params_changed)
        self.ref_width.trace_add("write", self._on_optical_params_changed)
        # core thickness only affects the 3D vertical stack geometry (not the
        # 2D n_eff), so it re-syncs the auto stack directly.
        self.core_thickness.trace_add("write", lambda *_a: self._resync_auto_stack())

    def _on_optical_params_changed(self, *_args):
        self._update_index_from_sellmeier()

    def _safe_float(self, var: tk.DoubleVar, default: Optional[float] = None) -> Optional[float]:
        try:
            return float(var.get())
        except Exception:
            return default

    def _update_index_from_sellmeier(self):
        if self._updating_index_vars:
            return
        lam = self._safe_float(self.wavelength)
        delta_pct = self._safe_float(self.delta_percent)
        if lam is None or delta_pct is None or lam <= 0:
            return
        self._updating_index_vars = True
        try:
            n_clad = BPM1DSolver.silica_sellmeier_n(lam)
            n_core = BPM1DSolver.n_core_from_delta_percent(n_clad, delta_pct)
            neff = BPM1DSolver.slab_te0_neff(lam, n_core, n_clad, self._estimate_reference_width_um())
            self.n_clad.set(round(n_clad, 9))
            self.n_core.set(round(n_core, 9))
            self.n_eff.set(round(neff, 9))
        finally:
            self._updating_index_vars = False
        # Keep the 3D vertical stack in step with the optical params: if the
        # stack is still the auto-generated default, refresh its layer indices
        # so a mid-session wavelength/Δ change actually reaches the 3D solver.
        self._resync_auto_stack()

    # ------------------------------------------------------------------
    # UI builder
    # ------------------------------------------------------------------

    def _build_ui(self):
        container = ttk.Frame(self)
        container.pack(fill=tk.BOTH, expand=True)

        left_outer = ttk.Frame(container)
        left_outer.pack(side=tk.LEFT, fill=tk.Y, padx=(8, 4), pady=8)
        right = ttk.Frame(container)
        right.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=(4, 8), pady=8)

        self.left_canvas = tk.Canvas(left_outer, highlightthickness=0, width=390)
        self.left_scrollbar = ttk.Scrollbar(left_outer, orient=tk.VERTICAL, command=self.left_canvas.yview)
        self.left_canvas.configure(yscrollcommand=self.left_scrollbar.set)
        self.left_canvas.pack(side=tk.LEFT, fill=tk.Y, expand=False)
        self.left_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        left = ttk.Frame(self.left_canvas)
        self.left_canvas_window = self.left_canvas.create_window((0, 0), window=left, anchor="nw")
        left.bind("<Configure>", self._on_left_frame_configure)
        self.left_canvas.bind("<Configure>", self._on_left_canvas_configure)
        self.left_canvas.bind_all("<MouseWheel>", self._on_mousewheel, add="+")
        self.left_canvas.bind_all("<Button-4>", self._on_mousewheel_linux, add="+")
        self.left_canvas.bind_all("<Button-5>", self._on_mousewheel_linux, add="+")

        self._section_file(left)
        self._section_bpm(left)
        self._section_stack(left)
        self._section_sweep(left)
        self._section_ports(left)
        self._section_actions(left)

        self.fig = Figure(figsize=(15.0, 8.3), dpi=100)
        gs = self.fig.add_gridspec(2, 3, hspace=0.38, wspace=0.32,
                                   left=0.06, right=0.97, top=0.94, bottom=0.07)
        self.ax_layout = self.fig.add_subplot(gs[0, 0])
        self.ax_index  = self.fig.add_subplot(gs[1, 0])
        self.ax_prop   = self.fig.add_subplot(gs[:, 1])   # propagation: full height
        self.ax_out    = self.fig.add_subplot(gs[0, 2])   # output power
        self.ax_phase  = self.fig.add_subplot(gs[1, 2])   # output phase

        self.canvas = FigureCanvasTkAgg(self.fig, master=right)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self._clear_axes()

    def _on_left_frame_configure(self, _event=None):
        if hasattr(self, "left_canvas"):
            self.left_canvas.configure(scrollregion=self.left_canvas.bbox("all"))

    def _on_left_canvas_configure(self, event):
        if hasattr(self, "left_canvas_window"):
            self.left_canvas.itemconfigure(self.left_canvas_window, width=event.width)

    def _left_canvas_has_focus(self) -> bool:
        if not hasattr(self, "left_canvas"):
            return False
        widget = self.winfo_containing(self.winfo_pointerx(), self.winfo_pointery())
        while widget is not None:
            if widget == self.left_canvas:
                return True
            widget = widget.master
        return False

    def _on_mousewheel(self, event):
        if self._left_canvas_has_focus():
            delta = event.delta
            if delta != 0:
                self.left_canvas.yview_scroll(int(-delta / 120), "units")

    def _on_mousewheel_linux(self, event):
        if self._left_canvas_has_focus():
            if getattr(event, "num", None) == 4:
                self.left_canvas.yview_scroll(-3, "units")
            elif getattr(event, "num", None) == 5:
                self.left_canvas.yview_scroll(3, "units")

    # ------------------------------------------------------------------
    # Panel sections
    # ------------------------------------------------------------------

    def _section_file(self, parent):
        frm = ttk.LabelFrame(parent, text="GDS")
        frm.pack(fill=tk.X, pady=4)
        frm.columnconfigure(0, weight=1)
        frm.columnconfigure(1, weight=1)
        frm.columnconfigure(2, weight=0)
        frm.columnconfigure(3, weight=0)

        ttk.Entry(frm, textvariable=self.gds_path, width=42).grid(row=0, column=0, columnspan=4, padx=4, pady=(4, 2), sticky="ew")
        self.btn_open = ttk.Button(frm, text="Open GDS", command=self.open_gds)
        self.btn_open.grid(row=1, column=0, columnspan=4, padx=4, pady=(0, 4), sticky="ew")

        ttk.Label(frm, text="Cell(optional)").grid(row=2, column=0, sticky="w", padx=4)
        ttk.Entry(frm, textvariable=self.cell_name, width=12).grid(row=2, column=1, sticky="ew", padx=4)
        ttk.Label(frm, text="Layer").grid(row=2, column=2, sticky="w", padx=4)
        ttk.Entry(frm, textvariable=self.layer, width=8).grid(row=2, column=3, padx=4, sticky="ew")

        ttk.Label(frm, text="Datatype").grid(row=3, column=0, sticky="w", padx=4)
        ttk.Entry(frm, textvariable=self.datatype, width=8).grid(row=3, column=1, padx=4, sticky="ew")
        ttk.Label(frm, text=f"gdstk: {'OK' if self.loader.available else 'missing'}").grid(row=3, column=2, columnspan=2, sticky="w", padx=4)

    def _section_bpm(self, parent):
        frm = ttk.LabelFrame(parent, text="BPM Parameters")
        frm.pack(fill=tk.X, pady=4)
        ttk.Label(frm, text="Wavelength [um]").grid(row=0, column=0, sticky="w", padx=4, pady=2)
        ttk.Entry(frm, textvariable=self.wavelength, width=12).grid(row=0, column=1, sticky="ew", padx=4, pady=2)
        ttk.Label(frm, text="Index contrast [%]").grid(row=1, column=0, sticky="w", padx=4, pady=2)
        ttk.Entry(frm, textvariable=self.delta_percent, width=12).grid(row=1, column=1, sticky="ew", padx=4, pady=2)
        ttk.Label(frm, text="n_clad (Sellmeier)").grid(row=2, column=0, sticky="w", padx=4, pady=2)
        ttk.Entry(frm, textvariable=self.n_clad, width=12, state="readonly").grid(row=2, column=1, sticky="ew", padx=4, pady=2)
        ttk.Label(frm, text="n_core (material)").grid(row=3, column=0, sticky="w", padx=4, pady=2)
        ttk.Entry(frm, textvariable=self.n_core, width=12, state="readonly").grid(row=3, column=1, sticky="ew", padx=4, pady=2)
        ttk.Label(frm, text="n_eff (auto)").grid(row=4, column=0, sticky="w", padx=4, pady=2)
        ttk.Entry(frm, textvariable=self.n_eff, width=12, state="readonly").grid(row=4, column=1, sticky="ew", padx=4, pady=2)
        for label, var, row in [
            ("dx [um]", self.dx, 5),
            ("dz [um]", self.dz, 6),
            ("x margin [um]", self.x_margin, 7),
            ("save every N steps", self.save_every, 8),
        ]:
            ttk.Label(frm, text=label).grid(row=row, column=0, sticky="w", padx=4, pady=2)
            ttk.Entry(frm, textvariable=var, width=12).grid(row=row, column=1, sticky="ew", padx=4, pady=2)

    def _section_stack(self, parent):
        frm = ttk.LabelFrame(parent, text="3D  /  Vertical Layer Stack")
        frm.pack(fill=tk.X, pady=4)
        frm.columnconfigure(1, weight=1)
        frm.columnconfigure(3, weight=1)

        # Simulation mode selector
        mode_row = ttk.Frame(frm)
        mode_row.grid(row=0, column=0, columnspan=4, sticky="ew", padx=4, pady=(4, 2))
        ttk.Label(mode_row, text="Simulation mode:").pack(side=tk.LEFT)
        ttk.Radiobutton(mode_row, text="2D (EIM)", value="2D",
                        variable=self.sim_mode, command=self._on_mode_change).pack(side=tk.LEFT, padx=6)
        ttk.Radiobutton(mode_row, text="3D (x,y,z)", value="3D",
                        variable=self.sim_mode, command=self._on_mode_change).pack(side=tk.LEFT)

        # Polarization + input-field selectors (3D only)
        pol_row = ttk.Frame(frm)
        pol_row.grid(row=1, column=0, columnspan=4, sticky="ew", padx=4, pady=(2, 1))
        ttk.Label(pol_row, text="Polarization:").pack(side=tk.LEFT)
        for txt, val in [("Scalar", "scalar"), ("TE", "TE"), ("TM", "TM"), ("Full", "full")]:
            ttk.Radiobutton(pol_row, text=txt, value=val, variable=self.pol_mode).pack(side=tk.LEFT, padx=2)

        inp_row = ttk.Frame(frm)
        inp_row.grid(row=2, column=0, columnspan=4, sticky="ew", padx=4, pady=(1, 4))
        ttk.Label(inp_row, text="Input:").pack(side=tk.LEFT)
        ttk.Radiobutton(inp_row, text="Gaussian", value="Gaussian",
                        variable=self.input_mode).pack(side=tk.LEFT, padx=2)
        ttk.Radiobutton(inp_row, text="Eigenmode", value="Eigenmode",
                        variable=self.input_mode).pack(side=tk.LEFT, padx=2)
        ttk.Label(inp_row, text="  Launch(full):").pack(side=tk.LEFT)
        ttk.Radiobutton(inp_row, text="Ex", value="Ex",
                        variable=self.launch_comp).pack(side=tk.LEFT, padx=2)
        ttk.Radiobutton(inp_row, text="Ey", value="Ey",
                        variable=self.launch_comp).pack(side=tk.LEFT, padx=2)

        ov_row = ttk.Frame(frm)
        ov_row.grid(row=7, column=0, columnspan=4, sticky="ew", padx=4, pady=(0, 2))
        ttk.Checkbutton(ov_row, text="Output eigenmode-overlap loss (extra mode solve)",
                        variable=self.out_mode_overlap).pack(side=tk.LEFT)

        # 3D numerical params
        grid_params = [
            ("dy [um]", self.dy, 0), ("y margin [um]", self.y_margin, 1),
            ("PML width [um]", self.pml_width, 2), ("bg index", self.bg_index, 3),
            ("bend R [um] (0=str)", self.bend_radius, 4), ("bend ctr [um]", self.bend_center, 5),
        ]
        prow = ttk.Frame(frm)
        prow.grid(row=3, column=0, columnspan=4, sticky="ew", padx=2)
        for label, var, i in grid_params:
            ttk.Label(prow, text=label).grid(row=i // 2, column=(i % 2) * 2, sticky="w", padx=4, pady=1)
            ttk.Entry(prow, textvariable=var, width=9).grid(row=i // 2, column=(i % 2) * 2 + 1, sticky="ew", padx=4, pady=1)

        yrow = ttk.Frame(frm)
        yrow.grid(row=4, column=0, columnspan=4, sticky="ew", padx=2)
        for label, var, i in [("Input y center [um]", self.input_y_center, 0),
                              ("Input y waist [um]", self.input_y_waist, 1)]:
            ttk.Label(yrow, text=label).grid(row=0, column=i * 2, sticky="w", padx=4, pady=1)
            ttk.Entry(yrow, textvariable=var, width=9).grid(row=0, column=i * 2 + 1, sticky="ew", padx=4, pady=1)

        # Stack table
        cols = ("name", "gds", "ybot", "thick", "n", "patt")
        self.stack_tree = ttk.Treeview(frm, columns=cols, show="headings", height=4)
        headers = [("name", "Layer", 70), ("gds", "GDS L/D", 56), ("ybot", "y0", 44),
                   ("thick", "h", 40), ("n", "index", 60), ("patt", "patt", 40)]
        for key, txt, w in headers:
            self.stack_tree.heading(key, text=txt)
            self.stack_tree.column(key, width=w, anchor="center")
        self.stack_tree.grid(row=5, column=0, columnspan=4, sticky="ew", padx=4, pady=(4, 2))
        self.stack_tree.bind("<Double-1>", lambda e: self._stack_edit())

        btns = ttk.Frame(frm)
        btns.grid(row=6, column=0, columnspan=4, sticky="ew", padx=4, pady=(0, 4))
        ttk.Button(btns, text="Auto", width=7, command=self._stack_autofill).pack(side=tk.LEFT, expand=True, fill=tk.X)
        ttk.Button(btns, text="Add", width=7, command=self._stack_add).pack(side=tk.LEFT, expand=True, fill=tk.X, padx=2)
        ttk.Button(btns, text="Edit", width=7, command=self._stack_edit).pack(side=tk.LEFT, expand=True, fill=tk.X)
        ttk.Button(btns, text="Remove", width=7, command=self._stack_remove).pack(side=tk.LEFT, expand=True, fill=tk.X, padx=2)

        self._on_mode_change()

    # ------------------------------------------------------------------
    # Vertical-stack editor
    # ------------------------------------------------------------------

    def _on_mode_change(self):
        is3d = self.sim_mode.get() == "3D"
        # The Run button label reflects the active mode
        if hasattr(self, "btn_run"):
            self.btn_run.configure(text="Run 3D BPM" if is3d else "Run BPM")
        if is3d and not self.stack_layers:
            self._stack_autofill()

    def _refresh_stack_tree(self):
        if not hasattr(self, "stack_tree"):
            return
        self.stack_tree.delete(*self.stack_tree.get_children())
        for i, l in enumerate(self.stack_layers):
            gds = f"{l.gds_layer}/{l.gds_datatype}" if l.patterned and l.gds_layer is not None else "-"
            self.stack_tree.insert("", "end", iid=str(i), values=(
                l.name, gds, f"{l.y_bottom_um:.2f}", f"{l.thickness_um:.2f}",
                f"{l.n_index:.4f}", "Y" if l.patterned else "n"))

    def _stack_autofill(self):
        """Build a default substrate / core / top-clad stack from current params."""
        n_clad = self._safe_float(self.n_clad) or 1.444
        n_core = self._safe_float(self.n_core) or 1.4549
        h = self._safe_float(self.core_thickness) or 6.0
        layer = int(self.layer.get()) if str(self.layer.get()).strip() else 1
        dtype = int(self.datatype.get()) if str(self.datatype.get()).strip() else 0
        self.bg_index.set(round(n_clad, 6))
        stack = default_silica_stack(n_core, n_clad, core_thickness_um=h,
                                     gds_layer=layer, gds_datatype=dtype)
        self.stack_layers = list(stack.layers)
        self._stack_is_auto = True
        # default input beam centred on the core
        self.input_y_center.set(round(stack.core_center_um(), 4))
        self.input_y_waist.set(round(max(h / 2.2, self.dy.get()), 4))
        self._refresh_stack_tree()

    def _resync_auto_stack(self):
        """Rebuild the auto-generated stack from the current optical params
        (n_core / n_clad) and core thickness, so mid-session changes to those
        actually reach the 3D solver.  No-op once the user has manually edited
        the stack (so hand-set layers are never clobbered).  The input-beam
        placement is left untouched — only `_stack_autofill` seeds it."""
        if not self._stack_is_auto or not self.stack_layers:
            return
        n_clad = self._safe_float(self.n_clad)
        n_core = self._safe_float(self.n_core)
        h = self._safe_float(self.core_thickness)
        if n_clad is None or n_core is None or not h or h <= 0:
            return
        layer = int(self.layer.get()) if str(self.layer.get()).strip() else 1
        dtype = int(self.datatype.get()) if str(self.datatype.get()).strip() else 0
        self.bg_index.set(round(n_clad, 6))
        stack = default_silica_stack(n_core, n_clad, core_thickness_um=h,
                                     gds_layer=layer, gds_datatype=dtype)
        self.stack_layers = list(stack.layers)
        self._refresh_stack_tree()

    def _selected_stack_index(self) -> Optional[int]:
        sel = self.stack_tree.selection()
        if not sel:
            return None
        try:
            return int(sel[0])
        except Exception:
            return None

    def _stack_remove(self):
        idx = self._selected_stack_index()
        if idx is None or not (0 <= idx < len(self.stack_layers)):
            return
        del self.stack_layers[idx]
        self._stack_is_auto = False
        self._refresh_stack_tree()

    def _stack_add(self):
        self._stack_dialog(None)

    def _stack_edit(self):
        idx = self._selected_stack_index()
        if idx is None or not (0 <= idx < len(self.stack_layers)):
            return
        self._stack_dialog(idx)

    def _stack_dialog(self, idx: Optional[int]):
        """Popup to add/edit one StackLayer."""
        existing = self.stack_layers[idx] if idx is not None else None
        win = tk.Toplevel(self)
        win.title("Edit stack layer" if existing else "Add stack layer")
        win.resizable(False, False)
        win.transient(self)

        v_name = tk.StringVar(value=existing.name if existing else "layer")
        v_n = tk.StringVar(value=str(existing.n_index if existing else self._safe_float(self.n_core) or 1.45))
        v_thick = tk.StringVar(value=str(existing.thickness_um if existing else 6.0))
        v_ybot = tk.StringVar(value=str(existing.y_bottom_um if existing else 0.0))
        v_patt = tk.BooleanVar(value=existing.patterned if existing else True)
        v_layer = tk.StringVar(value=str(existing.gds_layer if (existing and existing.gds_layer is not None) else self.layer.get()))
        v_dtype = tk.StringVar(value=str(existing.gds_datatype if existing else self.datatype.get()))

        rows = [("Name", v_name), ("Index n", v_n), ("Thickness [um]", v_thick),
                ("y bottom [um]", v_ybot), ("GDS layer", v_layer), ("GDS datatype", v_dtype)]
        for r, (label, var) in enumerate(rows):
            ttk.Label(win, text=label).grid(row=r, column=0, sticky="w", padx=8, pady=3)
            ttk.Entry(win, textvariable=var, width=16).grid(row=r, column=1, sticky="ew", padx=8, pady=3)
        ttk.Checkbutton(win, text="Patterned (use GDS layer mask)", variable=v_patt).grid(
            row=len(rows), column=0, columnspan=2, sticky="w", padx=8, pady=4)

        def save():
            try:
                layer_val = int(float(v_layer.get())) if v_patt.get() else None
                new = StackLayer(
                    name=v_name.get().strip() or "layer",
                    n_index=float(v_n.get()),
                    thickness_um=float(v_thick.get()),
                    y_bottom_um=float(v_ybot.get()),
                    patterned=bool(v_patt.get()),
                    gds_layer=layer_val,
                    gds_datatype=int(float(v_dtype.get())) if v_dtype.get().strip() else 0,
                )
            except ValueError as e:
                messagebox.showerror("Invalid value", str(e), parent=win)
                return
            if existing is not None:
                self.stack_layers[idx] = new
            else:
                self.stack_layers.append(new)
            self._stack_is_auto = False
            self._refresh_stack_tree()
            win.destroy()

        btn = ttk.Frame(win)
        btn.grid(row=len(rows) + 1, column=0, columnspan=2, pady=8)
        ttk.Button(btn, text="OK", command=save).pack(side=tk.LEFT, padx=6)
        ttk.Button(btn, text="Cancel", command=win.destroy).pack(side=tk.LEFT, padx=6)

    def _section_sweep(self, parent):
        frm = ttk.LabelFrame(parent, text="Wavelength Sweep")
        frm.pack(fill=tk.X, pady=4)
        for r, (label, var) in enumerate([
            ("Start [um]", self.lambda_start),
            ("Stop [um]", self.lambda_stop),
            ("Step [um]", self.lambda_step),
        ]):
            ttk.Label(frm, text=label).grid(row=r, column=0, sticky="w", padx=4, pady=2)
            ttk.Entry(frm, textvariable=var, width=12).grid(row=r, column=1, sticky="ew", padx=4, pady=2)
        ttk.Checkbutton(frm, text="Material dispersion (silica, 3D sweep)",
                        variable=self.sweep_dispersion).grid(
            row=3, column=0, columnspan=2, sticky="w", padx=4, pady=(4, 0))
        self.btn_sweep = ttk.Button(frm, text="Sweep Monitor Loss", command=self.run_wavelength_sweep)
        self.btn_sweep.grid(row=4, column=0, columnspan=2, sticky="ew", padx=4, pady=(6, 4))

    def _section_ports(self, parent):
        frm = ttk.LabelFrame(parent, text="Input / Output")
        frm.pack(fill=tk.X, pady=4)
        for r, (label, var) in enumerate([
            ("Input center [um]", self.input_center),
            ("Input waist [um]", self.input_waist),
            ("Input angle [deg]", self.input_angle),
            ("Output center [um]", self.out_center),
            ("Output waist [um]", self.out_waist),
            ("Monitor width [um]", self.monitor_width),
        ]):
            ttk.Label(frm, text=label).grid(row=r, column=0, sticky="w", padx=4, pady=2)
            ttk.Entry(frm, textvariable=var, width=12).grid(row=r, column=1, sticky="ew", padx=4, pady=2)

        ttk.Label(frm, text="Detected input port").grid(row=6, column=0, sticky="w", padx=4, pady=(8, 2))
        self.input_port_combo = ttk.Combobox(frm, textvariable=self.input_port_choice, state="readonly", width=24, values=["auto"])
        self.input_port_combo.grid(row=6, column=1, sticky="ew", padx=4, pady=(8, 2))
        self.input_port_combo.bind("<<ComboboxSelected>>", self.on_select_input_port)

        ttk.Label(frm, text="Detected output port").grid(row=7, column=0, sticky="w", padx=4, pady=2)
        self.output_port_combo = ttk.Combobox(frm, textvariable=self.output_port_choice, state="readonly", width=24, values=["auto"])
        self.output_port_combo.grid(row=7, column=1, sticky="ew", padx=4, pady=2)
        self.output_port_combo.bind("<<ComboboxSelected>>", self.on_select_output_port)

        ttk.Checkbutton(frm, text="Auto-apply first detected ports on load", variable=self.auto_apply_ports).grid(
            row=8, column=0, columnspan=2, sticky="w", padx=4, pady=(6, 2)
        )

        btn_row = ttk.Frame(frm)
        btn_row.grid(row=9, column=0, columnspan=2, sticky="ew", padx=4, pady=(8, 2))
        self.btn_detect = ttk.Button(btn_row, text="Auto-detect Ports", command=self.detect_ports)
        self.btn_detect.pack(side=tk.LEFT, fill=tk.X, expand=True)
        self.btn_apply_ports = ttk.Button(btn_row, text="Apply Auto Ports", command=self.apply_auto_ports)
        self.btn_apply_ports.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(6, 0))

    def _section_actions(self, parent):
        frm = ttk.LabelFrame(parent, text="Actions")
        frm.pack(fill=tk.X, pady=4)
        self.btn_build = ttk.Button(frm, text="Build Index Grid", command=self.build_index)
        self.btn_build.pack(fill=tk.X, padx=4, pady=3)
        self.btn_run = ttk.Button(frm, text="Run BPM", command=self.run_bpm)
        self.btn_run.pack(fill=tk.X, padx=4, pady=3)
        self.btn_save = ttk.Button(frm, text="Save Result JSON", command=self.save_json)
        self.btn_save.pack(fill=tk.X, padx=4, pady=3)
        self.progress_bar = ttk.Progressbar(frm, variable=self.progress_value, maximum=100.0, mode="determinate")
        self.progress_bar.pack(fill=tk.X, padx=4, pady=(8, 3))
        ttk.Label(frm, textvariable=self.progress_text, wraplength=340, justify=tk.LEFT).pack(fill=tk.X, padx=4, pady=(0, 6))
        ttk.Label(frm, textvariable=self.result_text, wraplength=340, justify=tk.LEFT).pack(fill=tk.X, padx=4, pady=6)

    # ------------------------------------------------------------------
    # Axes helpers
    # ------------------------------------------------------------------

    def _clear_axes(self):
        # Remove accumulated colorbars first
        for attr in ("_cbar_index", "_cbar_prop"):
            cbar = getattr(self, attr, None)
            if cbar is not None:
                try:
                    cbar.remove()
                except Exception:
                    pass
                setattr(self, attr, None)
        for cbar in getattr(self, "_cbars_3d", []):
            try:
                cbar.remove()
            except Exception:
                pass
        self._cbars_3d = []
        for ax, title in [
            (self.ax_layout, "Layout"),
            (self.ax_index,  "Index map"),
            (self.ax_prop,   "Propagation |ψ|²"),
            (self.ax_out,    "Output power"),
            (self.ax_phase,  "Output phase"),
        ]:
            ax.clear()
            ax.grid(True, alpha=0.25)
            ax.set_title(title)
        self.canvas.draw_idle()

    # ------------------------------------------------------------------
    # Progress / busy state
    # ------------------------------------------------------------------

    def set_progress(self, fraction: float, text: str):
        fraction = max(0.0, min(1.0, float(fraction)))
        self.progress_bar.configure(mode="determinate")
        self.progress_value.set(100.0 * fraction)
        self.progress_text.set(text)
        self.update_idletasks()

    def reset_progress(self, text: str = "Idle"):
        self.progress_bar.configure(mode="determinate")
        self.progress_value.set(0.0)
        self.progress_text.set(text)
        self.update_idletasks()

    def set_busy(self, busy: bool, text: Optional[str] = None, indeterminate: bool = False):
        self.busy = busy
        state = tk.DISABLED if busy else tk.NORMAL
        for btn in [self.btn_open, self.btn_detect, self.btn_apply_ports,
                    self.btn_build, self.btn_run, self.btn_save, self.btn_sweep]:
            btn.configure(state=state)
        if busy:
            if indeterminate:
                self.progress_bar.configure(mode="indeterminate")
                self.progress_bar.start(10)
            else:
                self.progress_bar.configure(mode="determinate")
                self.progress_bar.stop()
            if text:
                self.progress_text.set(text)
        else:
            self.progress_bar.stop()
            self.progress_bar.configure(mode="determinate")
            if text:
                self.progress_text.set(text)
        self.update_idletasks()

    # ------------------------------------------------------------------
    # Async worker
    # ------------------------------------------------------------------

    def _poll_worker_queue(self):
        try:
            while True:
                kind, payload = self.worker_queue.get_nowait()
                if kind == "progress":
                    frac, text = payload
                    self.set_progress(frac, text)
                elif kind == "done":
                    self.set_busy(False)
                    callback, result = payload
                    callback(result)
                elif kind == "error":
                    self.set_busy(False)
                    title, msg = payload
                    self.reset_progress(title)
                    messagebox.showerror(title, msg)
        except queue.Empty:
            pass
        self.after(80, self._poll_worker_queue)

    def _run_async(self, work_fn, done_fn, busy_text: str, indeterminate: bool = False):
        if self.busy:
            return
        self.set_busy(True, busy_text, indeterminate=indeterminate)

        def progress(frac: float, text: str):
            self.worker_queue.put(("progress", (frac, text)))

        def runner():
            try:
                result = work_fn(progress)
                self.worker_queue.put(("done", (done_fn, result)))
            except Exception as e:
                self.worker_queue.put(("error", ("Task error", str(e))))

        self.worker_thread = threading.Thread(target=runner, daemon=True)
        self.worker_thread.start()

    # ------------------------------------------------------------------
    # Solver helpers
    # ------------------------------------------------------------------

    def _selected_port_width(self, side: str) -> Optional[float]:
        ports = self.detected_input_ports if side == "input" else self.detected_output_ports
        choice = self.input_port_choice.get() if side == "input" else self.output_port_choice.get()
        prefix = "in" if side == "input" else "out"
        if choice.startswith(prefix):
            try:
                idx = int(choice.split(":", 1)[0][len(prefix):]) - 1
                if 0 <= idx < len(ports):
                    return float(ports[idx].width_um)
            except Exception:
                pass
        if ports:
            return float(ports[0].width_um)
        return None

    def _estimate_reference_width_um(self) -> float:
        width = self._selected_port_width("input")
        if width is None or width <= 0:
            width = max(2.0 * float(self.input_waist.get()), float(self.dx.get()))
        return float(width)

    def _snapshot_solver_params(self, wavelength_um: Optional[float] = None) -> dict:
        lam = float(self.wavelength.get() if wavelength_um is None else wavelength_um)
        delta_pct = float(self.delta_percent.get())
        dx = float(self.dx.get())
        dz = float(self.dz.get())
        x_margin = float(self.x_margin.get())
        ref_width = float(self._estimate_reference_width_um())
        n_clad = BPM1DSolver.silica_sellmeier_n(lam)
        n_core = BPM1DSolver.n_core_from_delta_percent(n_clad, delta_pct)
        return {
            "wavelength_um": lam,
            "delta_percent": delta_pct,
            "dx_um": dx,
            "dz_um": dz,
            "x_margin_um": x_margin,
            "reference_width_um": ref_width,
            "n_clad": n_clad,
            "n_core": n_core,
        }

    @staticmethod
    def _make_solver_from_params(params: dict) -> BPM1DSolver:
        return BPM1DSolver(
            wavelength_um=float(params["wavelength_um"]),
            n_core=float(params["n_core"]),
            n_clad=float(params["n_clad"]),
            dx_um=float(params["dx_um"]),
            dz_um=float(params["dz_um"]),
            x_margin_um=float(params["x_margin_um"]),
            reference_width_um=float(params["reference_width_um"]),
        )

    def current_solver(self, wavelength_um: Optional[float] = None) -> BPM1DSolver:
        params = self._snapshot_solver_params(wavelength_um)
        solver = self._make_solver_from_params(params)
        self.n_eff.set(round(solver.n_eff, 9))
        return solver

    def get_reference_width_um(self) -> Optional[float]:
        widths = []
        try:
            val = float(self.ref_width.get())
            if val > 0:
                widths.append(val)
        except Exception:
            pass
        choice = self.input_port_choice.get()
        if choice.startswith("in"):
            idx = int(choice.split(":", 1)[0][2:]) - 1
            if 0 <= idx < len(self.detected_input_ports):
                widths.append(float(self.detected_input_ports[idx].width_um))
        elif self.detected_input_ports:
            widths.append(float(self.detected_input_ports[0].width_um))
        return float(np.mean(widths)) if widths else None

    def _sync_neff_display_from_solver(self, solver: BPM1DSolver):
        self._updating_index_vars = True
        try:
            self.n_clad.set(round(solver.n_clad, 9))
            self.n_core.set(round(solver.n_core, 9))
            self.n_eff.set(round(solver.n_eff, 9))
        finally:
            self._updating_index_vars = False

    def _ensure_core_mask(
        self,
        solver: Optional[BPM1DSolver] = None,
        progress_cb: Optional[Callable[[float, str], None]] = None,
    ):
        if self.geometry_data is None:
            raise RuntimeError("Open a GDS file first.")
        if self.core_mask is not None and self.x is not None and self.z is not None:
            return self.x, self.z, self.core_mask
        if solver is None:
            solver = self.current_solver()
        x, z, core_mask = solver.build_core_mask_grid(self.geometry_data, progress_cb=progress_cb)
        self.x, self.z, self.core_mask = x, z, core_mask
        return x, z, core_mask

    def _build_n_map_for_current_wavelength(self):
        if self.core_mask is None:
            return None
        solver = self.current_solver()
        self.n_map, _ = solver.build_index_from_core_mask(self.x, self.core_mask)
        self.n_eff.set(round(solver.n_eff, 9))
        return self.n_map

    # ------------------------------------------------------------------
    # GDS open
    # ------------------------------------------------------------------

    def open_gds(self):
        # Always open the file dialog so the user can choose/change the file
        init_dir = os.path.dirname(self.gds_path.get().strip()) or os.path.expanduser("~")
        path = filedialog.askopenfilename(
            filetypes=[("GDS files", "*.gds *.GDS"), ("All files", "*.*")],
            initialdir=init_dir,
        )
        if not path:
            return
        self.gds_path.set(path)

        cell_name = self.cell_name.get().strip() or None
        layer = int(self.layer.get())
        datatype = int(self.datatype.get())
        solver_params = self._snapshot_solver_params()

        def work(progress_cb):
            progress_cb(0.05, "Scanning GDS layers/specs...")
            scan = self.loader.scan(path, cell_name=cell_name)
            progress_cb(0.25, "Loading selected layer/datatype...")
            chosen = (layer, datatype)
            try:
                geom = self.loader.load(path, layer=layer, datatype=datatype, cell_name=cell_name)
            except Exception:
                spec_pairs = [(ly, dt) for ly, dt, _ in scan.specs]
                if len(spec_pairs) == 1:
                    chosen = spec_pairs[0]
                    geom = self.loader.load(path, layer=chosen[0], datatype=chosen[1], cell_name=cell_name)
                else:
                    raise
            progress_cb(0.55, "Auto-detecting ports...")
            solver = self._make_solver_from_params(solver_params)
            left_ports, right_ports = solver.detect_ports(geom)
            progress_cb(1.0, "GDS load complete")
            return {
                "scan": scan,
                "geom": geom,
                "chosen": chosen,
                "left_ports": sorted(left_ports, key=lambda p: p.center_x_um, reverse=True),
                "right_ports": sorted(right_ports, key=lambda p: p.center_x_um, reverse=True),
            }

        def done(result):
            self.last_scan_info = result["scan"]
            chosen_layer, chosen_dtype = result["chosen"]
            self.layer.set(chosen_layer)
            self.datatype.set(chosen_dtype)
            self.geometry_data = result["geom"]
            self.detected_input_ports = result["left_ports"]
            self.detected_output_ports = result["right_ports"]
            self.core_mask = None
            self.x = None
            self.z = None
            self.n_map = None
            self.psi_out = None
            self.total_power = None
            self.power_map = None
            self.z_samples = None
            self.last_sweep = None
            self.result3d = None
            self._update_port_combos()
            if self.detected_input_ports:
                self.ref_width.set(round(float(self.detected_input_ports[0].width_um), 3))
            if self.auto_apply_ports.get():
                self.apply_auto_ports()
            self._sync_neff_display_from_solver(self.current_solver())
            self._plot_layout()
            for ax, title in [(self.ax_index, "Index map"), (self.ax_prop, "Propagation |ψ|²"),
                               (self.ax_out, "Output power"), (self.ax_phase, "Output phase")]:
                ax.clear(); ax.set_title(title); ax.grid(True, alpha=0.25)
            self.canvas.draw_idle()
            spec_text = ", ".join([f"L{ly}/D{dt} ({cnt})" for ly, dt, cnt in result["scan"].specs]) or "none"
            self.result_text.set(
                f"Loaded GDS polygons: {len(self.geometry_data.polygons)}\n"
                f"Available specs: {spec_text}\n"
                f"Detected input ports: {len(self.detected_input_ports)} | output ports: {len(self.detected_output_ports)}"
            )
            self.reset_progress("GDS loaded and ports auto-detected")

        self._run_async(work, done, "Opening GDS...", indeterminate=False)

    # ------------------------------------------------------------------
    # Plot helpers
    # ------------------------------------------------------------------

    def _plot_layout(self):
        self.ax_layout.clear()
        self.ax_layout.set_title("Layout (GDS X→z, Y→x)")
        self.ax_layout.set_xlabel("z [um]")
        self.ax_layout.set_ylabel("x [um]")
        if not self.geometry_data:
            self.canvas.draw_idle()
            return
        for poly in self.geometry_data.polygons:
            xy = np.asarray(poly.exterior.coords)
            patch = MplPolygon(xy, closed=True, fill=True, facecolor="steelblue", alpha=0.45, edgecolor="steelblue")
            self.ax_layout.add_patch(patch)
        minz, minx, maxz, maxx = self.geometry_data.bounds
        self.ax_layout.set_xlim(minz, maxz)
        self.ax_layout.set_ylim(minx - 2, maxx + 2)

        ic   = self.input_center.get()
        iw   = self.input_waist.get()
        oc   = self.out_center.get()
        ow   = self.out_waist.get()

        # Input: centre line + 1/e² waist band (입력면 전체에 걸쳐 반투명 띠)
        self.ax_layout.axhline(ic, linestyle="--", linewidth=1.2, color="tab:blue", label=f"Input ctr={ic:.3f}")
        self.ax_layout.axhspan(ic - iw, ic + iw, alpha=0.12, color="tab:blue",
                               label=f"Input waist ±{iw:.3f}")

        # Selected output: centre line + waist band
        self.ax_layout.axhline(oc, linestyle=":", linewidth=1.2, color="tab:orange", label=f"Out(sel) ctr={oc:.3f}")
        self.ax_layout.axhspan(oc - ow, oc + ow, alpha=0.12, color="tab:orange",
                               label=f"Out waist ±{ow:.3f}")

        # All detected output ports
        port_colors = ["tab:red", "tab:green", "tab:purple", "gold"]
        for i, p in enumerate(self.detected_output_ports):
            col = port_colors[i % len(port_colors)]
            self.ax_layout.axhline(p.center_x_um, color=col, linestyle=":",
                                   alpha=0.7, linewidth=1.2,
                                   label=f"out{i+1}: x={p.center_x_um:.3f}")

        self.ax_layout.legend(fontsize=7, loc="upper right")
        self.ax_layout.grid(True, alpha=0.25)
        self.canvas.draw_idle()

    def _plot_index(self):
        # Remove previous colorbar before clearing axes (prevents accumulation)
        if self._cbar_index is not None:
            try:
                self._cbar_index.remove()
            except Exception:
                pass
            self._cbar_index = None
        self.ax_index.clear()
        self.ax_index.set_title("Index map n(x, z)")
        self.ax_index.set_xlabel("z [um]")
        self.ax_index.set_ylabel("x [um]")
        if self.n_map is None:
            self.canvas.draw_idle()
            return
        im = self.ax_index.imshow(
            self.n_map.T, origin="lower", aspect="auto",
            extent=[self.z[0], self.z[-1], self.x[0], self.x[-1]],
            cmap="plasma",
        )
        self._cbar_index = self.fig.colorbar(im, ax=self.ax_index, pad=0.02, fraction=0.046, label="n_eff")
        self.ax_index.axhline(self.input_center.get(), linestyle="--", linewidth=1, color="cyan",  label="input")
        self.ax_index.axhline(self.out_center.get(),   linestyle=":",  linewidth=1, color="yellow", label="out(sel)")
        self.ax_index.legend(fontsize=7, loc="upper right")
        self.canvas.draw_idle()

    def _plot_propagation(self):
        # Remove previous colorbar before clearing axes (prevents accumulation)
        if self._cbar_prop is not None:
            try:
                self._cbar_prop.remove()
            except Exception:
                pass
            self._cbar_prop = None
        self.ax_prop.clear()
        self.ax_prop.set_title("Propagation |ψ(x,z)|²")
        self.ax_prop.set_xlabel("z [um]")
        self.ax_prop.set_ylabel("x [um]")
        im = self.ax_prop.imshow(
            self.power_map.T, origin="lower", aspect="auto",
            extent=[self.z_samples[0], self.z_samples[-1], self.x[0], self.x[-1]],
            cmap="inferno",
        )
        self._cbar_prop = self.fig.colorbar(im, ax=self.ax_prop, pad=0.02, fraction=0.03, label="|ψ|²")

        ic = self.input_center.get()
        iw = self.input_waist.get()
        oc = self.out_center.get()
        ow = self.out_waist.get()

        # Input beam: centre + waist band
        self.ax_prop.axhline(ic, linestyle="--", linewidth=1.2, color="cyan", label=f"In ctr={ic:.3f}")
        self.ax_prop.axhspan(ic - iw, ic + iw, alpha=0.15, color="cyan",
                             label=f"In waist ±{iw:.3f}")

        # Selected output: centre + waist band
        self.ax_prop.axhline(oc, linestyle=":", linewidth=1.2, color="tab:orange", label=f"Out ctr={oc:.3f}")
        self.ax_prop.axhspan(oc - ow, oc + ow, alpha=0.15, color="tab:orange",
                             label=f"Out waist ±{ow:.3f}")

        # All detected output ports
        port_colors = ["tab:red", "tab:green", "tab:purple", "gold"]
        for i, p in enumerate(self.detected_output_ports):
            col = port_colors[i % len(port_colors)]
            self.ax_prop.axhline(p.center_x_um, linestyle=":", linewidth=1.0,
                                 color=col, alpha=0.7, label=f"out{i+1}")

        self.ax_prop.legend(fontsize=7, loc="upper right")
        self.canvas.draw_idle()

    def _plot_output(self, mode_out: np.ndarray, port_modes: list = None,
                     port_powers: list = None, port_ils: list = None,
                     input_power: float = 1.0):
        port_colors = ["tab:red", "tab:green", "tab:purple", "gold"]
        pwr = np.abs(self.psi_out) ** 2

        # ── ax_out : Output power profile ──────────────────────────────
        self.ax_out.clear()
        self.ax_out.set_title("Output |ψ|²  — per-port")
        self.ax_out.set_xlabel("x [um]")
        self.ax_out.set_ylabel("|ψ|²")
        self.ax_out.plot(self.x, pwr, color="tab:blue", linewidth=1.5, label="|ψ_out|²", zorder=5)

        # Selected output waist 범위 (1/e² 기준)
        oc = float(self.out_center.get())
        ow = float(self.out_waist.get())
        self.ax_out.axvspan(oc - ow, oc + ow, alpha=0.08, color="tab:orange",
                            label=f"Out waist ±{ow:.3f}")
        self.ax_out.axvline(oc - ow, linestyle="-.", linewidth=0.8, color="tab:orange", alpha=0.6)
        self.ax_out.axvline(oc + ow, linestyle="-.", linewidth=0.8, color="tab:orange", alpha=0.6)

        # Shade each port window and annotate IL
        ports = self.detected_output_ports
        for i, p in enumerate(ports):
            col = port_colors[i % len(port_colors)]
            hw  = max(p.width_um * 0.75,
                      float(self.monitor_width.get()) / 2.0)
            win = np.abs(self.x - p.center_x_um) <= hw
            self.ax_out.fill_between(self.x, 0, pwr, where=win,
                                     alpha=0.25, color=col, label=f"port{i+1} window")
            # Target mode overlay
            if port_modes and i < len(port_modes):
                pm_pwr = np.abs(port_modes[i]) ** 2 * pwr.max() * 0.8  # scale to same range
                self.ax_out.plot(self.x, pm_pwr, linestyle="--",
                                 color=col, linewidth=1.0, alpha=0.8)
            # IL annotation
            il_txt = f"IL={port_ils[i]:.2f} dB" if port_ils and i < len(port_ils) else ""
            pw_txt = f"({port_powers[i]/max(input_power,1e-30)*100:.1f}%)" if port_powers and i < len(port_powers) else ""
            self.ax_out.axvline(p.center_x_um, linestyle=":", linewidth=1.0, color=col, alpha=0.7)
            ymax = pwr.max() if pwr.max() > 0 else 1.0
            self.ax_out.text(p.center_x_um, ymax * (0.85 - 0.18 * i),
                             f"out{i+1}\n{il_txt}\n{pw_txt}",
                             ha="center", va="top", fontsize=7,
                             color=col, bbox=dict(fc="white", ec=col, alpha=0.7, pad=2))
        self.ax_out.legend(fontsize=7, loc="upper left")
        self.ax_out.grid(True, alpha=0.25)

        # ── ax_phase : Output phase profile ────────────────────────────
        self.ax_phase.clear()
        self.ax_phase.set_title("Output phase ∠ψ(x)")
        self.ax_phase.set_xlabel("x [um]")
        self.ax_phase.set_ylabel("Phase [deg]")

        # Unwrap phase; mask regions where power is below 1% of peak
        threshold = pwr.max() * 0.01
        phase_deg = np.degrees(np.angle(self.psi_out))
        # Unwrap in place on the valid region
        valid = pwr >= threshold
        phase_plot = np.full_like(phase_deg, np.nan)
        if np.any(valid):
            idx_valid = np.where(valid)[0]
            phase_unwrapped = np.degrees(np.unwrap(np.angle(self.psi_out[valid])))
            # Centre phase at 0 (subtract mean)
            phase_unwrapped -= np.mean(phase_unwrapped)
            phase_plot[idx_valid] = phase_unwrapped

        self.ax_phase.plot(self.x, phase_plot, color="tab:cyan", linewidth=1.5, label="Phase (unwrapped)")

        # Selected output waist 범위 표시
        oc_p = float(self.out_center.get())
        ow_p = float(self.out_waist.get())
        self.ax_phase.axvspan(oc_p - ow_p, oc_p + ow_p, alpha=0.08, color="tab:orange",
                              label=f"Out waist ±{ow_p:.3f}")
        self.ax_phase.axvline(oc_p - ow_p, linestyle="-.", linewidth=0.8, color="tab:orange", alpha=0.6)
        self.ax_phase.axvline(oc_p + ow_p, linestyle="-.", linewidth=0.8, color="tab:orange", alpha=0.6)

        for i, p in enumerate(ports):
            col = port_colors[i % len(port_colors)]
            hw  = max(p.width_um * 0.75, float(self.monitor_width.get()) / 2.0)
            # Per-port mean phase
            win = (np.abs(self.x - p.center_x_um) <= hw) & valid
            if np.any(win):
                ph_mean = float(np.nanmean(phase_plot[win]))
                self.ax_phase.axhline(ph_mean, linestyle="--", color=col, linewidth=1.0, alpha=0.8)
                self.ax_phase.axvline(p.center_x_um, linestyle=":", color=col, linewidth=1.0, alpha=0.7)
                self.ax_phase.text(p.center_x_um, ph_mean,
                                   f" out{i+1}\n {ph_mean:.1f}°",
                                   fontsize=7, color=col, va="bottom")
        self.ax_phase.legend(fontsize=7)
        self.ax_phase.grid(True, alpha=0.25)

        self.canvas.draw_idle()

    # ------------------------------------------------------------------
    # Port management
    # ------------------------------------------------------------------

    def _update_port_combos(self):
        input_vals = ["auto"] + [f"in{i+1}: x={p.center_x_um:.2f} um, w={p.width_um:.2f} um" for i, p in enumerate(self.detected_input_ports)]
        output_vals = ["auto"] + [f"out{i+1}: x={p.center_x_um:.2f} um, w={p.width_um:.2f} um" for i, p in enumerate(self.detected_output_ports)]
        self.input_port_combo["values"] = input_vals
        self.output_port_combo["values"] = output_vals
        self.input_port_choice.set(input_vals[1] if len(input_vals) > 1 else "auto")
        self.output_port_choice.set(output_vals[1] if len(output_vals) > 1 else "auto")

    def apply_auto_ports(self):
        # Push the currently combo-selected port coordinates into the entry fields.
        # Input port
        if self.detected_input_ports:
            val_in = self.input_port_choice.get()
            try:
                idx_in = int(val_in.split(":", 1)[0][2:]) - 1
                if not (0 <= idx_in < len(self.detected_input_ports)):
                    idx_in = 0
            except Exception:
                idx_in = 0
            pin = self.detected_input_ports[idx_in]
            self.input_center.set(round(float(pin.center_x_um), 3))
            self.input_waist.set(round(max(self.dx.get() * 2.0, 0.45 * pin.width_um), 3))
            self.ref_width.set(round(float(pin.width_um), 3))
        # Output port
        if self.detected_output_ports:
            val_out = self.output_port_choice.get()
            try:
                idx_out = int(val_out.split(":", 1)[0][3:]) - 1
                if not (0 <= idx_out < len(self.detected_output_ports)):
                    idx_out = 0
            except Exception:
                idx_out = 0
            pout = self.detected_output_ports[idx_out]
            self.out_center.set(round(float(pout.center_x_um), 3))
            self.out_waist.set(round(max(self.dx.get() * 2.0, 0.45 * pout.width_um), 3))
            self.monitor_width.set(round(max(self.monitor_width.get(), 1.4 * pout.width_um), 3))
        self._sync_neff_display_from_solver(self.current_solver())
        self._plot_layout()

    def on_select_input_port(self, _event=None):
        # Combo selection only updates ref_width (for n_eff calculation).
        # input_center / input_waist are NOT overwritten — user controls those freely.
        val = self.input_port_choice.get()
        if not val.startswith("in"):
            return
        idx = int(val.split(":", 1)[0][2:]) - 1
        if 0 <= idx < len(self.detected_input_ports):
            p = self.detected_input_ports[idx]
            self.ref_width.set(round(float(p.width_um), 3))
            self._sync_neff_display_from_solver(self.current_solver())
            self._plot_layout()

    def on_select_output_port(self, _event=None):
        # Combo selection only highlights the port on the layout.
        # out_center / out_waist / monitor_width are NOT overwritten.
        # Use "Apply Auto Ports" to explicitly push a port's coordinates.
        self._plot_layout()

    def detect_ports(self, silent: bool = False):
        if not self.geometry_data:
            if not silent:
                messagebox.showwarning("No GDS", "Open a GDS file first.")
            return

        solver_params = self._snapshot_solver_params()

        def work(progress_cb):
            progress_cb(0.05, "Detecting ports...")
            solver = self._make_solver_from_params(solver_params)
            left_ports, right_ports = solver.detect_ports(self.geometry_data)
            progress_cb(1.0, f"Detected {len(left_ports)} input port(s), {len(right_ports)} output port(s)")
            return {
                "left_ports": sorted(left_ports, key=lambda p: p.center_x_um, reverse=True),
                "right_ports": sorted(right_ports, key=lambda p: p.center_x_um, reverse=True),
            }

        def done(result):
            self.detected_input_ports = result["left_ports"]
            self.detected_output_ports = result["right_ports"]
            self._update_port_combos()
            if self.auto_apply_ports.get():
                self.apply_auto_ports()
            self._sync_neff_display_from_solver(self.current_solver())
            self._plot_layout()
            self.result_text.set(
                "Detected ports\n"
                + "Input: " + (", ".join([f"{p.center_x_um:.2f} um" for p in self.detected_input_ports]) or "none")
                + "\nOutput: " + (", ".join([f"{p.center_x_um:.2f} um" for p in self.detected_output_ports]) or "none")
            )
            self.reset_progress("Port detection complete")
            if not silent:
                messagebox.showinfo("Port detection", self.result_text.get())

        self._run_async(work, done, "Detecting ports...", indeterminate=False)

    # ------------------------------------------------------------------
    # Build index grid
    # ------------------------------------------------------------------

    def build_index(self):
        if not self.geometry_data:
            messagebox.showwarning("No GDS", "Open a GDS file first.")
            return

        solver_params = self._snapshot_solver_params()

        def work(progress_cb):
            solver = self._make_solver_from_params(solver_params)
            x, z, core_mask = self._ensure_core_mask(solver=solver, progress_cb=progress_cb)
            n_map, neff_slices = solver.build_index_from_core_mask(x, core_mask)
            return {
                "x": x, "z": z, "core_mask": core_mask, "n_map": n_map,
                "n_clad": solver.n_clad, "n_core": solver.n_core, "n_eff": solver.n_eff,
                "neff_slices": neff_slices,
            }

        def done(result):
            self.x = result["x"]
            self.z = result["z"]
            self.core_mask = result["core_mask"]
            self.n_map = result["n_map"]
            self._plot_index()
            self.result_text.set(
                f"Index grid built: nx={len(self.x)}, nz={len(self.z)}\n"
                f"n_clad={result['n_clad']:.6f}, n_core={result['n_core']:.6f}, n_eff={result['n_eff']:.6f}\n"
                f"Guide region index used in BPM = {result['n_eff']:.6f}"
            )
            self.set_progress(1.0, f"Index grid complete: nx={len(self.x)}, nz={len(self.z)}")

        self._run_async(work, done, "Building index grid...", indeterminate=False)

    # ------------------------------------------------------------------
    # Monitor metrics  (fix #4 — single source of truth, used by run_bpm)
    # ------------------------------------------------------------------

    def _compute_monitor_metrics(
        self,
        solver: BPM1DSolver,
        x: np.ndarray,
        psi_out: np.ndarray,
        mode_out: np.ndarray,
        field0: np.ndarray,
    ) -> dict:
        # fix #5: call as instance method, not BPM1DSolver.overlap_power(solver, ...)
        channel_power = solver.overlap_power(x, psi_out, mode_out)
        input_power = float(np.trapz(np.abs(field0) ** 2, x))
        total_out = float(np.trapz(np.abs(psi_out) ** 2, x))
        width = float(self.monitor_width.get())
        window = np.abs(x - float(self.out_center.get())) <= width / 2
        window_power = float(np.trapz(np.abs(psi_out[window]) ** 2, x[window])) if np.any(window) else 0.0
        il_total = 10 * np.log10(max(input_power, 1e-30) / max(total_out, 1e-30))
        il_channel = 10 * np.log10(max(input_power, 1e-30) / max(channel_power, 1e-30))
        il_monitor = 10 * np.log10(max(input_power, 1e-30) / max(window_power, 1e-30))
        coupling_db = -10 * np.log10(max(channel_power, 1e-30))
        return {
            "input_power": input_power,
            "total_out": total_out,
            "window_power": window_power,
            "channel_power": channel_power,
            "il_total": il_total,
            "il_channel": il_channel,
            "il_monitor": il_monitor,
            "coupling_db": coupling_db,
        }

    # ------------------------------------------------------------------
    # Run BPM
    # ------------------------------------------------------------------

    def run_bpm(self):
        if self.sim_mode.get() == "3D":
            return self.run_bpm_3d()
        if self.geometry_data is None:
            messagebox.showwarning("No GDS", "Open a GDS file first.")
            return

        solver_params = self._snapshot_solver_params()
        input_center = float(self.input_center.get())
        input_waist  = float(self.input_waist.get())
        input_angle  = float(self.input_angle.get())
        out_center   = float(self.out_center.get())
        out_waist    = float(self.out_waist.get())
        monitor_width = float(self.monitor_width.get())
        save_every    = max(1, int(self.save_every.get()))
        # Snapshot ALL detected output ports for multi-port measurement
        out_ports_snap = list(self.detected_output_ports)

        def work(progress_cb):
            solver = self._make_solver_from_params(solver_params)
            x, z, core_mask = self._ensure_core_mask(
                solver=solver,
                progress_cb=lambda f, t: progress_cb(0.25 * f, t),
            )
            n_map, _ = solver.build_index_from_core_mask(x, core_mask)
            field0 = solver.make_input_field(x, x0_um=input_center, waist_um=input_waist, angle_deg=input_angle)
            psi_out, total_power, power_map, z_samples = solver.run(
                x, z, n_map, field0,
                save_every=save_every,
                progress_cb=lambda f, t: progress_cb(0.25 + 0.75 * f, t),
                progress_prefix="Running BPM",
            )
            input_power = float(np.trapz(np.abs(field0) ** 2, x))
            total_out   = float(np.trapz(np.abs(psi_out) ** 2, x))
            il_total    = 10 * np.log10(max(input_power, 1e-30) / max(total_out, 1e-30))

            # ── Single selected output port (overlap + window) ──────────────
            mode_out      = solver.gaussian_mode(x, center_um=out_center, waist_um=out_waist)
            channel_power = solver.overlap_power(x, psi_out, mode_out)
            win_sel       = np.abs(x - out_center) <= monitor_width / 2.0
            window_power  = (float(np.trapz(np.abs(psi_out[win_sel]) ** 2, x[win_sel]))
                             if np.any(win_sel) else 0.0)
            il_monitor  = 10 * np.log10(max(input_power, 1e-30) / max(window_power,  1e-30))
            il_channel  = 10 * np.log10(max(input_power, 1e-30) / max(channel_power, 1e-30))
            coupling_db = -10 * np.log10(max(channel_power, 1e-30))

            # ── Multi-port: measure power in EVERY detected output port ──────
            # Use 1.5× waveguide width as the per-port integration window so
            # adjacent ports don't overlap.
            port_powers = []
            port_modes  = []
            for p in out_ports_snap:
                hw = max(p.width_um * 0.75, monitor_width / 2.0)  # half-window
                win = np.abs(x - p.center_x_um) <= hw
                pw  = (float(np.trapz(np.abs(psi_out[win]) ** 2, x[win]))
                       if np.any(win) else 0.0)
                wst = max(solver_params["dx_um"] * 2, 0.45 * p.width_um)
                pm  = solver.gaussian_mode(x, center_um=p.center_x_um, waist_um=wst)
                port_powers.append(pw)
                port_modes.append(pm)
            total_port_power = sum(port_powers)
            il_all_ports = 10 * np.log10(max(input_power, 1e-30) / max(total_port_power, 1e-30))

            # Build per-port summary lines
            port_lines = []
            for i, (p, pw) in enumerate(zip(out_ports_snap, port_powers)):
                ratio = pw / max(input_power, 1e-30)
                il_p  = -10 * np.log10(max(ratio, 1e-30))
                port_lines.append(
                    f"  out{i+1} x={p.center_x_um:.2f} um: "
                    f"P={pw:.5f}  ({ratio*100:.1f}%)  IL={il_p:.3f} dB"
                )
            port_summary = "\n".join(port_lines) if port_lines else "  (no output ports detected)"

            port_ils = []
            for pw in port_powers:
                port_ils.append(10 * np.log10(max(input_power, 1e-30) / max(pw, 1e-30)))

            return {
                "x": x, "z": z, "core_mask": core_mask, "n_map": n_map,
                "psi_out": psi_out, "total_power": total_power,
                "power_map": power_map, "z_samples": z_samples,
                "mode_out": mode_out, "port_modes": port_modes,
                "port_powers": port_powers, "port_ils": port_ils,
                "input_power": input_power,
                "summary": (
                    f"BPM complete\n"
                    f"Input power        = {input_power:.6f}\n"
                    f"Total output power = {total_out:.6f}  IL={il_total:.3f} dB\n"
                    f"\n"
                    f"── Per output port ──\n"
                    f"{port_summary}\n"
                    f"  ALL ports sum: P={total_port_power:.5f}  IL={il_all_ports:.3f} dB\n"
                    f"\n"
                    f"── Selected port (x={out_center:.2f} um) ──\n"
                    f"Monitor window IL  = {il_monitor:.3f} dB\n"
                    f"Target-mode IL     = {il_channel:.3f} dB\n"
                    f"Mismatch loss      = {coupling_db:.3f} dB"
                ),
            }

        def done(result):
            self.result3d = None   # this is a 2D result
            self.x = result["x"]
            self.z = result["z"]
            self.core_mask = result["core_mask"]
            self.n_map = result["n_map"]
            self.psi_out = result["psi_out"]
            self.total_power = result["total_power"]
            self.power_map = result["power_map"]
            self.z_samples = result["z_samples"]
            self._plot_index()
            self._plot_propagation()
            self._plot_output(
                result["mode_out"],
                port_modes=result.get("port_modes", []),
                port_powers=result.get("port_powers", []),
                port_ils=result.get("port_ils", []),
                input_power=result.get("input_power", 1.0),
            )
            self.result_text.set(result["summary"])
            self.set_progress(1.0, "BPM run complete")

        self._run_async(work, done, "Running BPM...", indeterminate=False)

    # ------------------------------------------------------------------
    # 3D BPM
    # ------------------------------------------------------------------

    def _estimate_n_ref_3d(self, lam: float, stack: VerticalStack,
                           ref_w_um: Optional[float] = None) -> float:
        """Two-step effective index (vertical then lateral) for the 3D paraxial
        reference index — keeps the paraxial approximation accurate.

        ``ref_w_um`` may be supplied (snapshotted on the main thread) so this can
        be called from a worker thread without touching tk variables."""
        n_clad = stack.background_index
        patt = [l for l in stack.layers if l.patterned] or list(stack.layers)
        if not patt:
            return n_clad
        core = max(patt, key=lambda l: l.n_index)
        n_v = BPM1DSolver.slab_te0_neff(lam, core.n_index, n_clad, core.thickness_um)
        ref_w = self._estimate_reference_width_um() if ref_w_um is None else ref_w_um
        n_ref = BPM1DSolver.slab_te0_neff(lam, max(n_v, n_clad + 1e-6), n_clad, ref_w)
        return float(max(n_ref, n_clad + 1e-6))

    def run_bpm_3d(self):
        if self.geometry_data is None:
            messagebox.showwarning("No GDS", "Open a GDS file first.")
            return
        if not self.stack_layers:
            self._stack_autofill()
        if not any(l.patterned for l in self.stack_layers):
            messagebox.showwarning("No core layer",
                                    "The vertical stack has no patterned (waveguide) layer.\n"
                                    "Add one mapped to the GDS core layer.")
            return

        # Snapshot everything from tk vars on the main thread
        lam = float(self.wavelength.get())
        stack = VerticalStack(layers=list(self.stack_layers),
                              background_index=float(self.bg_index.get()))
        dx = float(self.dx.get()); dz = float(self.dz.get())
        dy = float(self.dy.get()); y_margin = float(self.y_margin.get())
        pml = float(self.pml_width.get())
        bend_r = float(self.bend_radius.get()); bend_c = float(self.bend_center.get())
        x_margin = float(self.x_margin.get())
        in_xc = float(self.input_center.get()); in_yc = float(self.input_y_center.get())
        in_wx = float(self.input_waist.get());  in_wy = float(self.input_y_waist.get())
        out_xc = float(self.out_center.get());  out_wx = float(self.out_waist.get())
        monitor_w = float(self.monitor_width.get())
        save_every = max(1, int(self.save_every.get()))
        pol = self.pol_mode.get()
        is_full = pol == "full"
        launch_comp = self.launch_comp.get()
        use_eigenmode = self.input_mode.get() == "Eigenmode"
        out_mode_ov = self.out_mode_overlap.get()
        n_ref = self._estimate_n_ref_3d(lam, stack)
        path = self.gds_path.get().strip()
        cell_name = self.cell_name.get().strip() or None
        primary_key = (int(self.layer.get()), int(self.datatype.get()))
        out_ports_snap = list(self.detected_output_ports)

        def work(progress_cb):
            solver1d = BPM1DSolver(wavelength_um=lam, n_core=stack.index_max(),
                                   n_clad=stack.background_index, dx_um=dx, dz_um=dz,
                                   x_margin_um=x_margin)
            # Common (x, z) grid + primary lateral mask
            x, z, primary_mask = solver1d.build_core_mask_grid(
                self.geometry_data, progress_cb=lambda f, t: progress_cb(0.20 * f, t))

            # Lateral mask for every patterned GDS layer in the stack
            masks: Dict[Tuple[int, int], np.ndarray] = {}
            keys = stack.mask_keys()
            for ki, key in enumerate(keys):
                if key == primary_key:
                    masks[key] = primary_mask
                    continue
                try:
                    geom = self.loader.load(path, layer=key[0], datatype=key[1], cell_name=cell_name)
                    _, _, m = solver1d.build_core_mask_grid(geom, x_in=x, z_in=z)
                    masks[key] = m
                except Exception:
                    masks[key] = np.zeros_like(primary_mask)  # layer absent → no core
                progress_cb(0.20 + 0.10 * (ki + 1) / max(1, len(keys)),
                            f"Lateral mask for GDS {key[0]}/{key[1]}")

            # Vertical grid + index profile pieces
            y = stack.make_y_grid(dy, y_margin)
            nx, ny = len(x), len(y)
            bg_profile = stack.background_profile(y)
            bands = stack.patterned_bands(y)
            bg2d = np.broadcast_to(bg_profile[None, :], (nx, ny))

            def index_fn(iz):
                n_xy = bg2d.copy()
                for key, band, n_val in bands:
                    lat = masks.get(key)
                    if lat is None:
                        continue
                    row = lat[iz]
                    if row.any():
                        n_xy[row[:, None] & band[None, :]] = n_val
                return n_xy

            # Eigenmode solving (when requested) uses a TE/TM seed; for full-vec
            # the seed polarization follows the launched component.
            mode_pol = pol
            if is_full:
                mode_pol = "TE" if launch_comp == "Ex" else "TM"

            mode_neff = None
            if use_eigenmode:
                # Solve the fundamental guided mode of the INPUT cross-section,
                # then propagate it with n_ref set to its own n_eff (most accurate).
                seed = BPM3DSolver(wavelength_um=lam, n_ref=n_ref, dx_um=dx, dy_um=dy,
                                   dz_um=dz, pml_um=pml, pol=mode_pol,
                                   bend_radius_um=bend_r, bend_center_um=bend_c)
                n2_in = index_fn(0) ** 2
                field0, mode_neff = seed.solve_mode(
                    x, y, n2_in, in_xc, in_yc, in_wx, in_wy,
                    progress_cb=lambda f, t: progress_cb(0.30 + 0.15 * f, t))
                n_ref_prop = mode_neff
                run_start = 0.45
            else:
                n_ref_prop = n_ref
                run_start = 0.30

            solver3d = BPM3DSolver(wavelength_um=lam, n_ref=n_ref_prop, dx_um=dx, dy_um=dy,
                                   dz_um=dz, pml_um=pml, pol=pol,
                                   bend_radius_um=bend_r, bend_center_um=bend_c)
            if not use_eigenmode:
                field0 = solver3d.gaussian_2d(x, y, in_xc, in_yc, in_wx, in_wy)

            prog = lambda f, t: progress_cb(run_start + (1 - run_start) * f, t)
            if is_full:
                # Launch the prepared field into the chosen component; the other
                # starts empty and is populated by cross-polarization coupling.
                zeros = np.zeros_like(field0)
                ex0, ey0 = (field0, zeros) if launch_comp == "Ex" else (zeros, field0)
                res = solver3d.run_fullvec(x, y, z, index_fn, ex0, ey0,
                                           save_every=save_every, progress_cb=prog)
            else:
                res = solver3d.run(x, y, z, index_fn, field0, save_every=save_every,
                                   progress_cb=prog)

            # ── metrics ──
            input_power = res.input_power
            pwr_xy = np.abs(res.psi_out) ** 2
            if is_full and res.psi_out_minor is not None:
                pwr_xy = pwr_xy + np.abs(res.psi_out_minor) ** 2
            total_out = float(np.trapz(np.trapz(pwr_xy, y, axis=1), x))
            il_total = 10 * np.log10(max(input_power, 1e-30) / max(total_out, 1e-30))

            mode_out = solver3d.gaussian_2d(x, y, out_xc, in_yc, out_wx, in_wy)
            channel_power = solver3d.overlap_power(res.psi_out, mode_out, x, y)
            il_channel = 10 * np.log10(max(input_power, 1e-30) / max(channel_power, 1e-30))

            # Optional: per-output-port modal loss matrix.  Solve each port's
            # eigenmode on the OUTPUT cross-section and couple the propagated
            # field into it.  Full-vector splits into TE (Ex→TE-mode) and
            # TM (Ey→TM-mode) columns, exposing cross-polarization routing.
            il_db = lambda e: -10 * np.log10(max(e, 1e-30))
            mode_ov_line = ""
            modal_matrix = None
            if out_mode_ov and out_ports_snap:
                n2_out = index_fn(len(z) - 1) ** 2
                n2_clad = float(n2_out.min())
                rows = []
                # Isolate one port's guide so the mode solver returns THAT port's
                # mode (not the global fundamental of a multi-port cross-section).
                def _port_index(cx, p):
                    hw = max(p.width_um * 1.5, out_wx * 2.0, monitor_w)
                    near = np.abs(x - cx) <= hw
                    return np.where(near[:, None], n2_out, n2_clad)
                if is_full:
                    te_solver = BPM3DSolver(wavelength_um=lam, n_ref=n_ref_prop, dx_um=dx,
                                            dy_um=dy, dz_um=dz, pml_um=pml, pol="TE",
                                            bend_radius_um=bend_r, bend_center_um=bend_c)
                    tm_solver = BPM3DSolver(wavelength_um=lam, n_ref=n_ref_prop, dx_um=dx,
                                            dy_um=dy, dz_um=dz, pml_um=pml, pol="TM",
                                            bend_radius_um=bend_r, bend_center_um=bend_c)
                    for i, p in enumerate(out_ports_snap):
                        progress_cb(0.97, f"Modal matrix: port {i+1}/{len(out_ports_snap)}")
                        n2_p = _port_index(p.center_x_um, p)
                        m_te, nte = te_solver.solve_mode(x, y, n2_p, p.center_x_um, in_yc,
                                                         out_wx, in_wy, max_iter=500)
                        m_tm, ntm = tm_solver.solve_mode(x, y, n2_p, p.center_x_um, in_yc,
                                                         out_wx, in_wy, max_iter=500)
                        e_te = solver3d.overlap_power(res.psi_out, m_te, x, y) / max(input_power, 1e-30)
                        e_tm = solver3d.overlap_power(res.psi_out_minor, m_tm, x, y) / max(input_power, 1e-30)
                        rows.append({"port": i + 1, "x_um": p.center_x_um, "neff_TE": nte,
                                     "neff_TM": ntm, "eff_TE": e_te, "eff_TM": e_tm,
                                     "IL_TE_dB": il_db(e_te), "IL_TM_dB": il_db(e_tm)})
                    lines = ["── Output modal-loss matrix  (rows=ports, cols=pol) ──",
                             "  port    x[um]     TE: IL (eff)        TM: IL (eff)"]
                    for r in rows:
                        lines.append(f"  out{r['port']:<2d}  {r['x_um']:7.2f}   "
                                     f"{r['IL_TE_dB']:7.3f}dB ({r['eff_TE']*100:5.1f}%)   "
                                     f"{r['IL_TM_dB']:7.3f}dB ({r['eff_TM']*100:5.1f}%)")
                    mode_ov_line = "\n".join(lines) + "\n"
                else:
                    for i, p in enumerate(out_ports_snap):
                        progress_cb(0.97, f"Modal matrix: port {i+1}/{len(out_ports_snap)}")
                        n2_p = _port_index(p.center_x_um, p)
                        m_out, mneff = solver3d.solve_mode(x, y, n2_p, p.center_x_um, in_yc,
                                                          out_wx, in_wy, max_iter=500)
                        eff = solver3d.overlap_power(res.psi_out, m_out, x, y) / max(input_power, 1e-30)
                        rows.append({"port": i + 1, "x_um": p.center_x_um, "neff": mneff,
                                     "eff": eff, "IL_dB": il_db(eff)})
                    lines = ["── Output modal coupling (per port) ──"]
                    for r in rows:
                        lines.append(f"  out{r['port']:<2d} x={r['x_um']:7.2f}  n_eff={r['neff']:.5f}  "
                                     f"IL={r['IL_dB']:.3f} dB ({r['eff']*100:.1f}%)")
                    mode_ov_line = "\n".join(lines) + "\n"
                modal_matrix = rows

            # per output port: integrate total transverse power over the window
            port_lines = []
            for i, p in enumerate(out_ports_snap):
                hw = max(p.width_um * 0.75, monitor_w / 2.0)
                win = np.abs(x - p.center_x_um) <= hw
                pw = float(np.trapz(np.trapz(pwr_xy[win], y, axis=1), x[win])) if win.any() else 0.0
                ratio = pw / max(input_power, 1e-30)
                il_p = -10 * np.log10(max(ratio, 1e-30))
                port_lines.append(f"  out{i+1} x={p.center_x_um:.2f} um: "
                                  f"P={pw:.5f} ({ratio*100:.1f}%) IL={il_p:.3f} dB")
            port_summary = "\n".join(port_lines) if port_lines else "  (no output ports detected)"

            mode_line = (f"input = eigenmode, n_eff = {mode_neff:.6f}\n"
                         if mode_neff is not None else "input = Gaussian\n")
            pol_line = f"pol={pol}"
            if bend_r:
                pol_line += f", bend R={bend_r:.1f}um"
            if is_full:
                pol_line += (f", launch={launch_comp}, "
                             f"pol-conversion = {res.conversion*100:.3f}% (-> minor comp)")
            return {
                "res": res, "x": x, "y": y, "z": z, "n_ref": n_ref_prop,
                "mode_neff": mode_neff, "n_input": index_fn(0), "is_full": is_full,
                "modal_matrix": modal_matrix,
                "summary": (
                    f"3D BPM complete (FD-ADI + PML, {pol_line})\n"
                    f"grid nx={nx}, ny={ny}, nz={len(z)}   n_ref={n_ref_prop:.5f}\n"
                    f"{mode_line}"
                    f"Input power        = {input_power:.6f}\n"
                    f"Total output power = {total_out:.6f}  IL={il_total:.3f} dB\n"
                    f"\n── Per output port ──\n{port_summary}\n"
                    f"\n── Selected port (x={out_xc:.2f} um) ──\n"
                    f"Target-mode overlap IL = {il_channel:.3f} dB\n"
                    f"{mode_ov_line}"
                ),
            }

        def done(result):
            self.result3d = result["res"]
            self.mode_neff3d = result.get("mode_neff")
            self.modal_matrix3d = result.get("modal_matrix")
            self.x = result["x"]
            self.y = result["y"]
            self.z = result["z"]
            self.psi_out = result["res"].psi_out
            self.total_power = result["res"].total_power
            self._plot_3d_results(result["res"], result["n_input"])
            self.result_text.set(result["summary"])
            self.set_progress(1.0, "3D BPM run complete")

        self._run_async(work, done, "Running 3D BPM...", indeterminate=False)

    def _plot_3d_results(self, res, n_input):
        # Clear all colorbars/axes first to avoid accumulation
        self._clear_axes()
        self._plot_layout()

        cm = "plasma"
        # Input vertical index cross-section n(x, y)
        im0 = self.ax_index.imshow(n_input.T, origin="lower", aspect="auto",
                                   extent=[res.x[0], res.x[-1], res.y[0], res.y[-1]], cmap=cm)
        self.ax_index.set_title("Input index n(x, y)")
        self.ax_index.set_xlabel("x [um]"); self.ax_index.set_ylabel("y [um]")
        self._cbars_3d.append(self.fig.colorbar(im0, ax=self.ax_index, pad=0.02, fraction=0.046))

        # Top view |psi|^2 (z, x)
        im1 = self.ax_prop.imshow(res.topview.T, origin="lower", aspect="auto",
                                  extent=[res.z_samples[0], res.z_samples[-1], res.x[0], res.x[-1]],
                                  cmap="inferno")
        self.ax_prop.set_title("Top view  ∫|ψ|²dy  (z, x)")
        self.ax_prop.set_xlabel("z [um]"); self.ax_prop.set_ylabel("x [um]")
        for i, p in enumerate(self.detected_output_ports):
            self.ax_prop.axhline(p.center_x_um, linestyle=":", linewidth=1.0,
                                 color=["tab:red", "tab:green", "tab:purple", "gold"][i % 4],
                                 alpha=0.7, label=f"out{i+1}")
        if self.detected_output_ports:
            self.ax_prop.legend(fontsize=7, loc="upper right")
        self._cbars_3d.append(self.fig.colorbar(im1, ax=self.ax_prop, pad=0.02, fraction=0.03))

        full = getattr(res, "psi_out_minor", None) is not None
        ext_xy = [res.x[0], res.x[-1], res.y[0], res.y[-1]]

        if full:
            # Full-vector: show the converted minor component |Ey|² here instead
            # of the side view, so the polarization mixing is visible.
            ey = np.abs(res.psi_out_minor) ** 2
            im2 = self.ax_phase.imshow(ey.T, origin="lower", aspect="auto",
                                       extent=ext_xy, cmap="magma")
            conv = (res.conversion or 0.0) * 100
            self.ax_phase.set_title(f"Output minor |E_y(x,y)|²  (conv {conv:.2f}%)")
            self.ax_phase.set_xlabel("x [um]"); self.ax_phase.set_ylabel("y [um]")
        else:
            # Side view |psi|^2 (z, y)
            im2 = self.ax_phase.imshow(res.sideview.T, origin="lower", aspect="auto",
                                       extent=[res.z_samples[0], res.z_samples[-1], res.y[0], res.y[-1]],
                                       cmap="inferno")
            self.ax_phase.set_title("Side view  ∫|ψ|²dx  (z, y)")
            self.ax_phase.set_xlabel("z [um]"); self.ax_phase.set_ylabel("y [um]")
        self._cbars_3d.append(self.fig.colorbar(im2, ax=self.ax_phase, pad=0.02, fraction=0.03))

        # Output cross-section: |psi|^2 (scalar/TE/TM) or |Ex|^2 (full-vector)
        pwr = np.abs(res.psi_out) ** 2
        im3 = self.ax_out.imshow(pwr.T, origin="lower", aspect="auto",
                                 extent=ext_xy, cmap="viridis")
        self.ax_out.set_title("Output |E_x(x,y)|²" if full else "Output |ψ(x, y)|²")
        self.ax_out.set_xlabel("x [um]"); self.ax_out.set_ylabel("y [um]")
        self.ax_out.axvline(float(self.out_center.get()), linestyle=":", color="white", linewidth=1)
        self._cbars_3d.append(self.fig.colorbar(im3, ax=self.ax_out, pad=0.02, fraction=0.046))

        self.canvas.draw_idle()

    # ------------------------------------------------------------------
    # 3D wavelength sweep (IL + polarization conversion vs wavelength)
    # ------------------------------------------------------------------

    def run_wavelength_sweep_3d(self):
        if self.geometry_data is None:
            messagebox.showwarning("No GDS", "Open a GDS file first.")
            return
        if not self.stack_layers:
            self._stack_autofill()

        stack = VerticalStack(layers=list(self.stack_layers),
                              background_index=float(self.bg_index.get()))
        dx = float(self.dx.get()); dz = float(self.dz.get())
        dy = float(self.dy.get()); y_margin = float(self.y_margin.get())
        pml = float(self.pml_width.get()); x_margin = float(self.x_margin.get())
        in_xc = float(self.input_center.get()); in_yc = float(self.input_y_center.get())
        in_wx = float(self.input_waist.get());  in_wy = float(self.input_y_waist.get())
        out_xc = float(self.out_center.get())
        monitor_w = float(self.monitor_width.get())
        save_every = max(1, int(self.save_every.get()))
        pol = self.pol_mode.get()
        is_full = pol == "full"
        launch_comp = self.launch_comp.get()
        path = self.gds_path.get().strip()
        cell_name = self.cell_name.get().strip() or None
        primary_key = (int(self.layer.get()), int(self.datatype.get()))
        disp = self.sweep_dispersion.get()
        lam_ref = float(self.wavelength.get())
        ref_w = self._estimate_reference_width_um()   # snapshot on main thread
        bend_r = float(self.bend_radius.get()); bend_c = float(self.bend_center.get())

        def scale_stack(st, r):
            return VerticalStack(
                background_index=st.background_index * r,
                layers=[StackLayer(l.name, l.n_index * r, l.thickness_um, l.y_bottom_um,
                                   l.patterned, l.gds_layer, l.gds_datatype) for l in st.layers])

        def work(progress_cb):
            lambdas = self._sweep_wavelength_values()
            nlam = len(lambdas)
            # geometry (wavelength-independent) built once
            solver1d = BPM1DSolver(wavelength_um=float(lambdas[0]), n_core=stack.index_max(),
                                   n_clad=stack.background_index, dx_um=dx, dz_um=dz,
                                   x_margin_um=x_margin)
            x, z, primary_mask = solver1d.build_core_mask_grid(
                self.geometry_data, progress_cb=lambda f, t: progress_cb(0.10 * f, t))
            masks: Dict[Tuple[int, int], np.ndarray] = {}
            for key in stack.mask_keys():
                if key == primary_key:
                    masks[key] = primary_mask
                    continue
                try:
                    geom = self.loader.load(path, layer=key[0], datatype=key[1], cell_name=cell_name)
                    _, _, m = solver1d.build_core_mask_grid(geom, x_in=x, z_in=z)
                    masks[key] = m
                except Exception:
                    masks[key] = np.zeros_like(primary_mask)
            y = stack.make_y_grid(dy, y_margin)   # y geometry is index-independent
            nx, ny = len(x), len(y)

            def make_index_fn(cur_stack):
                bands = cur_stack.patterned_bands(y)
                bg2d = np.broadcast_to(cur_stack.background_profile(y)[None, :], (nx, ny))

                def index_fn(iz):
                    n_xy = bg2d.copy()
                    for key, band, n_val in bands:
                        lat = masks.get(key)
                        if lat is None:
                            continue
                        row = lat[iz]
                        if row.any():
                            n_xy[row[:, None] & band[None, :]] = n_val
                    return n_xy
                return index_fn

            n_ref_lam0 = BPM1DSolver.silica_sellmeier_n(lam_ref)

            il_total, il_sel, conv = [], [], []
            for i, lam in enumerate(lambdas):
                # silica material dispersion: scale all stack indices by the
                # Sellmeier ratio n_silica(λ)/n_silica(λ_ref) (uniform approximation)
                if disp:
                    ratio = BPM1DSolver.silica_sellmeier_n(float(lam)) / n_ref_lam0
                    cur_stack = scale_stack(stack, ratio)
                else:
                    cur_stack = stack
                index_fn = make_index_fn(cur_stack)
                n_ref = self._estimate_n_ref_3d(float(lam), cur_stack, ref_w)
                solver3d = BPM3DSolver(wavelength_um=float(lam), n_ref=n_ref, dx_um=dx,
                                       dy_um=dy, dz_um=dz, pml_um=pml, pol=pol,
                                       bend_radius_um=bend_r, bend_center_um=bend_c)
                base = 0.10 + 0.90 * (i / max(nlam, 1))
                span = 0.90 / max(nlam, 1)
                cb = lambda f, t, b=base, s=span: progress_cb(b + s * f, f"λ={lam:.4f}: {t}")
                field0 = solver3d.gaussian_2d(x, y, in_xc, in_yc, in_wx, in_wy)
                if is_full:
                    zeros = np.zeros_like(field0)
                    ex0, ey0 = (field0, zeros) if launch_comp == "Ex" else (zeros, field0)
                    res = solver3d.run_fullvec(x, y, z, index_fn, ex0, ey0,
                                               save_every=save_every, progress_cb=cb)
                    pwr_xy = np.abs(res.psi_out) ** 2 + np.abs(res.psi_out_minor) ** 2
                    conv.append(res.conversion * 100.0)
                else:
                    res = solver3d.run(x, y, z, index_fn, field0,
                                       save_every=save_every, progress_cb=cb)
                    pwr_xy = np.abs(res.psi_out) ** 2
                    conv.append(0.0)
                inp = res.input_power
                tot = float(np.trapz(np.trapz(pwr_xy, y, axis=1), x))
                il_total.append(10 * np.log10(max(inp, 1e-30) / max(tot, 1e-30)))
                win = np.abs(x - out_xc) <= monitor_w / 2.0
                psel = float(np.trapz(np.trapz(pwr_xy[win], y, axis=1), x[win])) if win.any() else 0.0
                il_sel.append(10 * np.log10(max(inp, 1e-30) / max(psel, 1e-30)))

            return {
                "lambdas": lambdas, "il_total": np.array(il_total),
                "il_sel": np.array(il_sel), "conversion": np.array(conv),
                "is_full": is_full, "pol": pol, "dispersion": disp,
            }

        def done(result):
            self.last_sweep3d = result
            self._show_sweep3d_plot(result)
            self.result_text.set(
                f"3D sweep complete [{result['pol']}]: {len(result['lambdas'])} wavelengths\n"
                f"Total IL {result['il_total'].min():.2f}–{result['il_total'].max():.2f} dB"
                + (f"\nPol-conversion {result['conversion'].min():.2f}–{result['conversion'].max():.2f}%"
                   if result['is_full'] else "")
            )
            self.set_progress(1.0, "3D wavelength sweep complete")

        self._run_async(work, done, "Running 3D wavelength sweep...", indeterminate=False)

    def _show_sweep3d_plot(self, result):
        if self.sweep_window is None or not self.sweep_window.winfo_exists():
            self.sweep_window = tk.Toplevel(self)
            self.sweep_window.title("3D Sweep — Loss / Conversion vs Wavelength")
            self.sweep_window.geometry("820x540")
            fig = Figure(figsize=(7.9, 5.0), dpi=100)
            ax = fig.add_subplot(111)
            canvas = FigureCanvasTkAgg(fig, master=self.sweep_window)
            canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
            self.sweep_window.fig = fig
            self.sweep_window.ax = ax
            self.sweep_window.canvas = canvas
        ax = self.sweep_window.ax
        ax.clear()
        lam = result["lambdas"]
        ax.plot(lam, result["il_total"], "-o", ms=3, label="Total insertion loss [dB]")
        ax.plot(lam, result["il_sel"], "-s", ms=3, label="Selected-port loss [dB]")
        ax.set_xlabel("Wavelength [um]"); ax.set_ylabel("Loss [dB]")
        ax.grid(True, alpha=0.25)
        disp_txt = "silica dispersion" if result.get("dispersion") else "fixed index"
        title = f"3D sweep [{result['pol']}], Gaussian input, {disp_txt}"
        if result["is_full"]:
            ax2 = ax.twinx()
            ax2.plot(lam, result["conversion"], "-^", color="tab:red", ms=3,
                     label="Pol-conversion [%]")
            ax2.set_ylabel("Polarization conversion [%]", color="tab:red")
            ax2.tick_params(axis="y", labelcolor="tab:red")
            ax2.legend(loc="upper right", fontsize=8)
        ax.set_title(title)
        ax.legend(loc="upper left", fontsize=8)
        self.sweep_window.canvas.draw_idle()
        self.sweep_window.lift()

    # ------------------------------------------------------------------
    # Wavelength sweep
    # ------------------------------------------------------------------

    def _sweep_wavelength_values(self) -> np.ndarray:
        start = float(self.lambda_start.get())
        stop = float(self.lambda_stop.get())
        step = float(self.lambda_step.get())
        if step <= 0:
            raise ValueError("Wavelength step must be > 0.")
        if stop < start:
            raise ValueError("Wavelength stop must be >= start.")
        count = int(np.floor((stop - start) / step + 1e-12)) + 1
        vals = start + step * np.arange(count)
        if vals[-1] < stop - 1e-12:
            vals = np.append(vals, stop)
        return np.round(vals, 12)

    def run_wavelength_sweep(self):
        if self.sim_mode.get() == "3D":
            return self.run_wavelength_sweep_3d()
        if self.geometry_data is None:
            messagebox.showwarning("No GDS", "Open a GDS file first.")
            return

        base_solver_params = self._snapshot_solver_params()
        input_center = float(self.input_center.get())
        input_waist = float(self.input_waist.get())
        input_angle = float(self.input_angle.get())
        out_center = float(self.out_center.get())
        out_waist = float(self.out_waist.get())
        monitor_width = float(self.monitor_width.get())
        out_ports_snap = list(self.detected_output_ports)  # snapshot for multi-port sweep

        def work(progress_cb):
            base_solver = self._make_solver_from_params(base_solver_params)
            x, z, core_mask = self._ensure_core_mask(
                solver=base_solver,
                progress_cb=lambda f, t: progress_cb(0.15 * f, t),
            )
            lambdas = self._sweep_wavelength_values()
            monitor_losses, channel_losses, total_losses = [], [], []
            nclads, ncores, neffs = [], [], []
            nlam = len(lambdas)

            for i, lam in enumerate(lambdas):
                sp = dict(base_solver_params)
                sp["wavelength_um"] = float(lam)
                sp["n_clad"] = BPM1DSolver.silica_sellmeier_n(float(lam))
                sp["n_core"] = BPM1DSolver.n_core_from_delta_percent(sp["n_clad"], sp["delta_percent"])
                solver = self._make_solver_from_params(sp)
                n_map, _ = solver.build_index_from_core_mask(x, core_mask)
                field0 = solver.make_input_field(x, x0_um=input_center, waist_um=input_waist, angle_deg=input_angle)
                base = 0.15 + 0.85 * (i / max(nlam, 1))
                span = 0.85 / max(nlam, 1)
                psi_out, _, _, _ = solver.run(
                    x, z, n_map, field0,
                    save_every=max(999999, len(z) + 1),
                    progress_cb=lambda f, t, b=base, s=span, lv=lam: progress_cb(b + s * f, f"Sweep {lv:.4f} um | {t}"),
                    progress_prefix=f"Sweep {lam:.4f} um",
                )
                mode_out = solver.gaussian_mode(x, center_um=out_center, waist_um=out_waist)
                channel_power = solver.overlap_power(x, psi_out, mode_out)
                input_power = float(np.trapz(np.abs(field0) ** 2, x))
                total_out = float(np.trapz(np.abs(psi_out) ** 2, x))
                window = np.abs(x - out_center) <= monitor_width / 2.0
                window_power = float(np.trapz(np.abs(psi_out[window]) ** 2, x[window])) if np.any(window) else 0.0
                # Multi-port: sum power in ALL detected output ports
                all_port_pw = 0.0
                for p in out_ports_snap:
                    hw = max(p.width_um * 0.75, monitor_width / 2.0)
                    win = np.abs(x - p.center_x_um) <= hw
                    all_port_pw += float(np.trapz(np.abs(psi_out[win]) ** 2, x[win])) if np.any(win) else 0.0
                # Use all-port power as the primary monitor metric for splitters
                eff_monitor_pw = all_port_pw if all_port_pw > 0 else window_power

                monitor_losses.append(10 * np.log10(max(input_power, 1e-30) / max(eff_monitor_pw, 1e-30)))
                channel_losses.append(10 * np.log10(max(input_power, 1e-30) / max(channel_power, 1e-30)))
                total_losses.append(10 * np.log10(max(input_power, 1e-30) / max(total_out, 1e-30)))
                nclads.append(solver.n_clad)
                ncores.append(solver.n_core)
                neffs.append(solver.n_eff)   # fix #3: only one append per iteration

            return {
                "lambdas": np.array(lambdas, dtype=float),
                "monitor_losses": np.array(monitor_losses, dtype=float),
                "channel_losses": np.array(channel_losses, dtype=float),
                "total_losses": np.array(total_losses, dtype=float),
                "nclads": np.array(nclads, dtype=float),
                "ncores": np.array(ncores, dtype=float),
                "neffs": np.array(neffs, dtype=float),
                "x": x, "z": z, "core_mask": core_mask,
            }

        def done(result):
            self.last_sweep = result
            self.x = result["x"]
            self.z = result["z"]
            self.core_mask = result["core_mask"]
            self._show_sweep_plot(result)
            idx_best = int(np.argmin(result["monitor_losses"]))
            self.result_text.set(
                f"Sweep complete: {len(result['lambdas'])} points\n"
                f"Best monitor loss = {result['monitor_losses'][idx_best]:.3f} dB"
                f" at {result['lambdas'][idx_best]:.4f} um\n"
                f"n_clad={result['nclads'][idx_best]:.6f}, n_core={result['ncores'][idx_best]:.6f}"
            )
            self.set_progress(1.0, "Wavelength sweep complete")

        self._run_async(work, done, "Running wavelength sweep...", indeterminate=False)

    def _show_sweep_plot(self, result):
        if self.sweep_window is None or not self.sweep_window.winfo_exists():
            self.sweep_window = tk.Toplevel(self)
            self.sweep_window.title("Monitor Loss vs Wavelength")
            self.sweep_window.geometry("780x520")
            fig = Figure(figsize=(7.5, 4.8), dpi=100)
            ax = fig.add_subplot(111)
            canvas = FigureCanvasTkAgg(fig, master=self.sweep_window)
            canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
            self.sweep_window.fig = fig
            self.sweep_window.ax = ax
            self.sweep_window.canvas = canvas
        ax = self.sweep_window.ax
        ax.clear()
        ax.plot(result["lambdas"], result["monitor_losses"], label="Monitor-port loss [dB]")
        ax.plot(result["lambdas"], result["channel_losses"], label="Target-channel loss [dB]")
        ax.plot(result["lambdas"], result["total_losses"], label="Total insertion loss [dB]")
        ax.set_xlabel("Wavelength [um]")
        ax.set_ylabel("Loss [dB]")
        ax.set_title("Loss spectrum")
        ax.grid(True, alpha=0.25)
        ax.legend()
        self.sweep_window.canvas.draw_idle()
        self.sweep_window.lift()

    # ------------------------------------------------------------------
    # Save JSON
    # ------------------------------------------------------------------

    def save_json(self):
        if self.psi_out is None and self.last_sweep is None and self.result3d is None:
            messagebox.showwarning("No result", "Run BPM or wavelength sweep first.")
            return
        path = filedialog.asksaveasfilename(defaultextension=".json", filetypes=[("JSON", "*.json")])
        if not path:
            return
        data = {
            "wavelength_um": float(self.wavelength.get()),
            "delta_percent": float(self.delta_percent.get()),
            "n_clad": float(self.n_clad.get()),
            "n_core": float(self.n_core.get()),
            "n_eff": float(self.n_eff.get()),
            "reference_width_um": float(self.ref_width.get()),
            "detected_input_ports": [p.__dict__ for p in self.detected_input_ports],
            "detected_output_ports": [p.__dict__ for p in self.detected_output_ports],
        }
        if self.result3d is not None:
            r = self.result3d
            data["mode"] = "3D"
            data["polarization"] = self.pol_mode.get()
            data["input_field"] = self.input_mode.get()
            data["bend_radius_um"] = float(self.bend_radius.get())
            data["bend_center_um"] = float(self.bend_center.get())
            if self.mode_neff3d is not None:
                data["eigenmode_neff"] = float(self.mode_neff3d)
            data["vertical_stack"] = [l.__dict__ for l in self.stack_layers]
            data.update({
                "x_um": r.x.tolist(),
                "y_um": r.y.tolist(),
                "z_samples_um": r.z_samples.tolist(),
                "output_power_xy": (np.abs(r.psi_out) ** 2).tolist(),
                "topview_zx": r.topview.tolist(),
                "sideview_zy": r.sideview.tolist(),
                "total_power_vs_z": r.total_power.tolist(),
            })
            if getattr(r, "psi_out_minor", None) is not None:
                data["launch_component"] = self.launch_comp.get()
                data["polarization_conversion"] = float(r.conversion or 0.0)
                data["output_power_minor_xy"] = (np.abs(r.psi_out_minor) ** 2).tolist()
            if self.modal_matrix3d:
                data["modal_loss_matrix"] = self.modal_matrix3d
        elif self.psi_out is not None:
            data["mode"] = "2D"
            data.update({
                "x_um": self.x.tolist(),
                "z_um": self.z.tolist(),
                "output_power": (np.abs(self.psi_out) ** 2).tolist(),
                "output_amplitude": np.abs(self.psi_out).tolist(),
                "output_phase_rad": np.angle(self.psi_out).tolist(),
                "total_power_vs_z": self.total_power.tolist(),
            })
        if self.last_sweep is not None:
            data["wavelength_sweep"] = {
                "lambda_um": self.last_sweep["lambdas"].tolist(),
                "monitor_loss_db": self.last_sweep["monitor_losses"].tolist(),
                "channel_loss_db": self.last_sweep["channel_losses"].tolist(),
                "total_loss_db": self.last_sweep["total_losses"].tolist(),
                "n_clad": self.last_sweep["nclads"].tolist(),
                "n_core": self.last_sweep["ncores"].tolist(),
                "n_eff": self.last_sweep["neffs"].tolist(),
            }
        with open(path, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2)
        self.result_text.set(f"Saved: {path}")


if __name__ == "__main__":
    app = BPMGui()
    app.mainloop()