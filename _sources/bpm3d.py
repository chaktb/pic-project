"""
3D scalar BPM solver + vertical layer-stack model
==================================================

Companion to ``bpm_gui_rev08.py``.  The existing ``BPM1DSolver`` propagates a
field on the (x, z) plane and collapses the vertical (y) direction with the
effective-index method — in photonics terms that is a *2D* simulation.

This module adds a genuine *3D* solver, ``BPM3DSolver``, that propagates a
scalar field ``psi(x, y, z)`` with two transverse dimensions:

  * x  — lateral (GDS Y axis)
  * y  — vertical / out-of-plane (defined by a layer stack, NOT in the GDS)
  * z  — propagation (GDS X axis)

Because a GDS file has no vertical information, the vertical refractive-index
profile is supplied by a :class:`VerticalStack`: each entry maps a *GDS layer*
to a vertical slab (y-bottom, thickness, index).  Where that GDS layer's
polygons exist laterally (a 2D core mask on the (z, x) grid, exactly what the
1D/2D tool already produces), the slab's index is painted into the y range;
elsewhere the background index fills in.

Numerics
--------
Paraxial (Fresnel) scalar equation, field convention ``E = psi * exp(+i*beta_ref*z)``
matching ``BPM1DSolver``::

    2 i beta_ref dpsi/dz = (d2/dx2 + d2/dy2) psi + (k0^2 n^2 - beta_ref^2) psi

advanced with the Peaceman–Rachford **ADI** scheme (one transverse direction
implicit per half-step, Crank–Nicolson, unconditionally stable, 2nd order in z).
The transverse Laplacian uses a **stretched-coordinate PML** so radiation leaves
the window without reflection.  Tridiagonal systems are solved with a batched
Thomas algorithm (pure numpy — no scipy dependency).

The solver only ever materialises the 2D index slice ``n(x, y)`` for the current
z step, so memory stays at O(nx*ny) for the field plus O(nz*nx) bool masks.
"""

from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Vertical stack model
# ---------------------------------------------------------------------------

@dataclass
class StackLayer:
    """One vertical slab of the layer stack.

    A layer occupies the vertical range ``[y_bottom_um, y_bottom_um + thickness_um]``.

    * If ``patterned`` is True the layer's index is applied only where the
      mapped GDS layer's polygons exist laterally (a 2D core mask); elsewhere
      the lower (background) layers show through.
    * If ``patterned`` is False the layer fills the whole (x, z) cross-section
      for its y range — use this for substrate / buffer / top cladding.
    """

    name: str
    n_index: float
    thickness_um: float
    y_bottom_um: float = 0.0
    patterned: bool = False
    gds_layer: Optional[int] = None
    gds_datatype: int = 0

    @property
    def y_top_um(self) -> float:
        return self.y_bottom_um + self.thickness_um

    def mask_key(self) -> Optional[Tuple[int, int]]:
        if not self.patterned or self.gds_layer is None:
            return None
        return (int(self.gds_layer), int(self.gds_datatype))


@dataclass
class VerticalStack:
    """Ordered set of :class:`StackLayer` plus a background fill index.

    ``background_index`` fills any y that no layer covers (typically the
    cladding).  Patterned layers are painted in list order, so a later entry
    overrides an earlier one where they overlap in (x, y, z).
    """

    layers: List[StackLayer] = field(default_factory=list)
    background_index: float = 1.444

    # -- y range / grid -------------------------------------------------
    def y_extent_um(self, margin_um: float = 0.0) -> Tuple[float, float]:
        if not self.layers:
            return (-1.0 - margin_um, 1.0 + margin_um)
        ymin = min(l.y_bottom_um for l in self.layers)
        ymax = max(l.y_top_um for l in self.layers)
        if ymax <= ymin:
            ymax = ymin + 1.0
        return (ymin - margin_um, ymax + margin_um)

    def make_y_grid(self, dy_um: float, margin_um: float) -> np.ndarray:
        ymin, ymax = self.y_extent_um(margin_um)
        # Symmetric grid centred on the stack centre so a centred mode stays centred.
        yc = 0.5 * (ymin + ymax)
        half = 0.5 * (ymax - ymin)
        n_half = max(1, int(np.ceil(half / dy_um)))
        return yc + dy_um * np.arange(-n_half, n_half + 1)

    def core_center_um(self) -> float:
        """Vertical centre of the highest-index patterned (or any) layer."""
        cand = [l for l in self.layers if l.patterned] or list(self.layers)
        if not cand:
            return 0.0
        core = max(cand, key=lambda l: l.n_index)
        return 0.5 * (core.y_bottom_um + core.y_top_um)

    def index_max(self) -> float:
        vals = [l.n_index for l in self.layers] + [self.background_index]
        return float(max(vals))

    def mask_keys(self) -> List[Tuple[int, int]]:
        keys = []
        for l in self.layers:
            k = l.mask_key()
            if k is not None and k not in keys:
                keys.append(k)
        return keys

    # -- index profile --------------------------------------------------
    def background_profile(self, y: np.ndarray) -> np.ndarray:
        """1D background index n0(y) (uniform layers only)."""
        n0 = np.full(y.shape, self.background_index, dtype=float)
        for l in self.layers:
            if l.patterned:
                continue
            sel = (y >= l.y_bottom_um) & (y < l.y_top_um)
            n0[sel] = l.n_index
        return n0

    def patterned_bands(self, y: np.ndarray) -> List[Tuple[Tuple[int, int], np.ndarray, float]]:
        """For each patterned layer return (mask_key, y-band boolean, n_index)."""
        bands = []
        for l in self.layers:
            k = l.mask_key()
            if k is None:
                continue
            band = (y >= l.y_bottom_um) & (y < l.y_top_um)
            bands.append((k, band, float(l.n_index)))
        return bands


def default_silica_stack(n_core: float, n_clad: float, core_thickness_um: float = 6.0,
                         gds_layer: int = 1, gds_datatype: int = 0) -> VerticalStack:
    """A simple buried symmetric core (one patterned GDS layer)."""
    return VerticalStack(
        background_index=n_clad,
        layers=[
            StackLayer("substrate", n_clad, thickness_um=8.0, y_bottom_um=-8.0,
                       patterned=False),
            StackLayer("core", n_core, thickness_um=core_thickness_um,
                       y_bottom_um=-0.5 * core_thickness_um, patterned=True,
                       gds_layer=gds_layer, gds_datatype=gds_datatype),
            StackLayer("top_clad", n_clad, thickness_um=8.0,
                       y_bottom_um=0.5 * core_thickness_um, patterned=False),
        ],
    )


# ---------------------------------------------------------------------------
# Batched tridiagonal (Thomas) solver
# ---------------------------------------------------------------------------

def thomas_batch(lower: np.ndarray, diag: np.ndarray, upper: np.ndarray,
                 rhs: np.ndarray) -> np.ndarray:
    """Solve tridiagonal systems along axis 0, vectorised over axis 1.

    lower/upper : (n,) or (n, m)   sub/super-diagonals (lower[0], upper[-1] unused)
    diag        : (n, m)           main diagonal (may vary per column)
    rhs         : (n, m)
    returns     : (n, m) solution
    """
    n = diag.shape[0]
    m = diag.shape[1]
    low = np.broadcast_to(lower.reshape(n, -1), (n, m))
    up = np.broadcast_to(upper.reshape(n, -1), (n, m))

    cp = np.empty((n, m), dtype=np.complex128)
    dp = np.empty((n, m), dtype=np.complex128)
    cp[0] = up[0] / diag[0]
    dp[0] = rhs[0] / diag[0]
    for i in range(1, n):
        denom = diag[i] - low[i] * cp[i - 1]
        cp[i] = up[i] / denom
        dp[i] = (rhs[i] - low[i] * dp[i - 1]) / denom
    x = np.empty((n, m), dtype=np.complex128)
    x[-1] = dp[-1]
    for i in range(n - 2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i + 1]
    return x


# ---------------------------------------------------------------------------
# 1D transverse operator with stretched-coordinate PML
# ---------------------------------------------------------------------------

def _pml_stretch(coord: np.ndarray, dh: float, pml_um: float,
                 sigma_max: float, m: float) -> Tuple[np.ndarray, np.ndarray]:
    """Return complex PML stretch factor s at nodes and at half-nodes (i+1/2).

    s = 1 + i * sigma_max * (d/L)^m  inside the PML layer of width ``pml_um``
    at each end, s = 1 in the interior.
    """
    lo = coord[0]
    hi = coord[-1]

    def stretch(c):
        s = np.ones_like(c, dtype=np.complex128)
        if pml_um <= 0:
            return s
        d_lo = (lo + pml_um) - c          # >0 inside the low-end PML
        d_hi = c - (hi - pml_um)          # >0 inside the high-end PML
        sel = d_lo > 0
        s[sel] += 1j * sigma_max * (d_lo[sel] / pml_um) ** m
        sel = d_hi > 0
        s[sel] += 1j * sigma_max * (d_hi[sel] / pml_um) ** m
        return s

    s_node = stretch(coord)
    s_half = stretch(coord + 0.5 * dh)    # s at i+1/2
    return s_node, s_half


def _laplacian_tridiag(coord: np.ndarray, dh: float, alpha: complex,
                       pml_um: float, sigma_max: float, m: float):
    """Tridiagonal coefficients of  alpha * d2/dh2  with stretched-coord PML.

    Operator (per line):  alpha * (1/s_i) d/dh[(1/s) d psi/dh]
    Returns (lower, diag, upper), each length len(coord), with Dirichlet ends.
    """
    n = len(coord)
    s_node, s_half_plus = _pml_stretch(coord, dh, pml_um, sigma_max, m)
    # s at i-1/2 is s_half_plus shifted by one (i.e. stretch evaluated at i-1/2)
    s_half_minus = np.empty_like(s_half_plus)
    s_half_minus[1:] = s_half_plus[:-1]
    s_half_minus[0] = s_half_plus[0]

    inv_dh2 = 1.0 / (dh * dh)
    lower = alpha * inv_dh2 / (s_node * s_half_minus)
    upper = alpha * inv_dh2 / (s_node * s_half_plus)
    diag = -alpha * inv_dh2 / s_node * (1.0 / s_half_plus + 1.0 / s_half_minus)

    # Dirichlet boundaries: no coupling outside the window.
    lower = lower.astype(np.complex128)
    upper = upper.astype(np.complex128)
    diag = diag.astype(np.complex128)
    lower[0] = 0.0
    upper[-1] = 0.0
    return lower, diag, upper


def _corrected_laplacian_tridiag(coord: np.ndarray, dh: float, alpha: complex,
                                 pml_um: float, sigma_max: float, m: float,
                                 n2: np.ndarray, axis: int):
    """Semi-vectorial transverse operator  alpha * d/dh[(1/n^2) d(n^2 psi)/dh]
    with stretched-coordinate PML, for the field component normal to interfaces
    along ``axis`` (TE -> axis=0/x, TM -> axis=1/y).

    Returns 2D (nx, ny) tridiagonal coefficients (lower, diag, upper) acting
    along ``axis``.  In a uniform region (n2 const) this reduces exactly to the
    scalar PML operator, so the PML layers (in cladding) are unaffected.

    Conservative flux form on a uniform grid:
        flux_{i+1/2} = (n2_{i+1} psi_{i+1} - n2_i psi_i) / (s_{i+1/2} * nbar2_{i+1/2})
        L psi_i = alpha/s_i * (flux_{i+1/2} - flux_{i-1/2}) / dh^2
    with nbar2 the arithmetic mean of neighbouring n^2.
    """
    n = len(coord)
    s_node, s_half_plus = _pml_stretch(coord, dh, pml_um, sigma_max, m)
    s_half_minus = np.empty_like(s_half_plus)
    s_half_minus[1:] = s_half_plus[:-1]
    s_half_minus[0] = s_half_plus[0]
    inv_dh2 = 1.0 / (dh * dh)

    # move the active axis to position 0 for a uniform treatment, then move back
    n2w = n2 if axis == 0 else n2.T          # (n, other)
    n2_plus = np.empty_like(n2w)             # n2 at i+1
    n2_minus = np.empty_like(n2w)            # n2 at i-1
    n2_plus[:-1] = n2w[1:]
    n2_plus[-1] = n2w[-1]
    n2_minus[1:] = n2w[:-1]
    n2_minus[0] = n2w[0]
    nbar_plus = 0.5 * (n2w + n2_plus)        # nbar2 at i+1/2
    nbar_minus = 0.5 * (n2w + n2_minus)      # nbar2 at i-1/2

    sN = s_node[:, None]
    sP = s_half_plus[:, None]
    sM = s_half_minus[:, None]

    lower = alpha * inv_dh2 / sN * (n2_minus / (sM * nbar_minus))
    upper = alpha * inv_dh2 / sN * (n2_plus / (sP * nbar_plus))
    diag = -alpha * inv_dh2 / sN * n2w * (1.0 / (sP * nbar_plus) + 1.0 / (sM * nbar_minus))

    lower = lower.astype(np.complex128)
    upper = upper.astype(np.complex128)
    diag = diag.astype(np.complex128)
    lower[0, :] = 0.0
    upper[-1, :] = 0.0
    if axis != 0:
        lower, diag, upper = lower.T, diag.T, upper.T
    return lower, diag, upper


def _as_axis_array(c: np.ndarray, axis: int) -> np.ndarray:
    """Broadcast a tridiagonal coefficient to 2D.

    Accepts either a 1D array along ``axis`` (line-constant operator) or an
    already-2D array (index-dependent, e.g. semi-vectorial TE/TM).
    """
    c = np.asarray(c)
    if c.ndim == 2:
        return c
    return c[:, None] if axis == 0 else c[None, :]


def _tridiag_matvec(lower: np.ndarray, diag: np.ndarray, upper: np.ndarray,
                    psi: np.ndarray, axis: int) -> np.ndarray:
    """Apply a tridiagonal operator along ``axis`` of psi.

    lower/diag/upper may be 1D (line-constant) or 2D (index-dependent).
    Dirichlet ends.
    """
    L = _as_axis_array(lower, axis)
    D = _as_axis_array(diag, axis)
    U = _as_axis_array(upper, axis)
    up_shift = np.zeros_like(psi)
    lo_shift = np.zeros_like(psi)
    if axis == 0:
        up_shift[:-1] = psi[1:]          # psi[i+1]
        lo_shift[1:] = psi[:-1]          # psi[i-1]
    else:
        up_shift[:, :-1] = psi[:, 1:]
        lo_shift[:, 1:] = psi[:, :-1]
    return L * lo_shift + D * psi + U * up_shift


def _thomas_axis(lower: np.ndarray, diag: np.ndarray, upper: np.ndarray,
                 rhs: np.ndarray, axis: int) -> np.ndarray:
    """Solve (lower, diag, upper) tridiagonal systems along ``axis`` of rhs.

    Coefficients may be 1D (along axis) or 2D.  diag and rhs are 2D.
    """
    if axis == 0:
        return thomas_batch(lower, diag, upper, rhs)
    lo = lower.T if np.ndim(lower) == 2 else lower
    up = upper.T if np.ndim(upper) == 2 else upper
    return thomas_batch(lo, diag.T, up, rhs.T).T


def _ddx(f: np.ndarray, dx: float) -> np.ndarray:
    """Central first derivative along x (axis 0), zero at the boundaries."""
    g = np.zeros_like(f)
    g[1:-1] = (f[2:] - f[:-2]) / (2.0 * dx)
    return g


def _ddy(f: np.ndarray, dy: float) -> np.ndarray:
    """Central first derivative along y (axis 1), zero at the boundaries."""
    g = np.zeros_like(f)
    g[:, 1:-1] = (f[:, 2:] - f[:, :-2]) / (2.0 * dy)
    return g


# ---------------------------------------------------------------------------
# 3D solver
# ---------------------------------------------------------------------------

@dataclass
class BPM3DResult:
    x: np.ndarray
    y: np.ndarray
    z_samples: np.ndarray
    psi_out: np.ndarray            # (nx, ny) complex final field (Ex for full-vec)
    total_power: np.ndarray        # (nz,) power vs z
    topview: np.ndarray            # (n_save, nx)   integral over y, vs z
    sideview: np.ndarray           # (n_save, ny)   integral over x, vs z
    input_power: float
    # full-vectorial extras (None for scalar/TE/TM)
    psi_out_minor: Optional[np.ndarray] = None   # Ey output component
    conversion: Optional[float] = None           # minor / total output power


class BPM3DSolver:
    def __init__(
        self,
        wavelength_um: float,
        n_ref: float,
        dx_um: float,
        dy_um: float,
        dz_um: float,
        pml_um: float = 2.0,
        pml_sigma: float = 4.0,
        pml_m: float = 3.0,
        pol: str = "scalar",
        bend_radius_um: float = 0.0,
        bend_center_um: float = 0.0,
    ):
        self.wavelength_um = float(wavelength_um)
        self.k0 = 2.0 * np.pi / self.wavelength_um
        self.n_ref = float(n_ref)
        self.beta_ref = self.k0 * self.n_ref
        self.dx_um = float(dx_um)
        self.dy_um = float(dy_um)
        self.dz_um = float(dz_um)
        self.pml_um = float(pml_um)
        self.pml_sigma = float(pml_sigma)
        self.pml_m = float(pml_m)
        # Bend: 0 (or None) means a straight guide.  A finite radius applies the
        # equivalent-straight-waveguide conformal index transform (see _bent).
        self.bend_radius_um = float(bend_radius_um or 0.0)
        self.bend_center_um = float(bend_center_um)
        pol = (pol or "scalar").upper()
        # "TE": dominant E along x (lateral) -> index correction along x
        # "TM": dominant E along y (vertical) -> index correction along y
        # "FULL": full-vectorial (Ex, Ey) with cross-polarization coupling
        self.pol = pol if pol in ("TE", "TM", "FULL") else "SCALAR"

    # -- grids ----------------------------------------------------------
    @staticmethod
    def make_symmetric_grid(lo: float, hi: float, dh: float) -> np.ndarray:
        c = 0.5 * (lo + hi)
        half = 0.5 * (hi - lo)
        n_half = max(1, int(np.ceil(half / dh)))
        return c + dh * np.arange(-n_half, n_half + 1)

    # -- bend (equivalent straight waveguide) ---------------------------
    def _bent(self, n2: np.ndarray, x: np.ndarray) -> np.ndarray:
        """Apply the constant-radius bend conformal transform to n².

        A bend of radius R is equivalent to a straight guide with index
        n_eq(x) = n(x)·exp((x-x0)/R)  ⇒  n_eq² = n²·exp(2(x-x0)/R).
        The lateral index ramp pushes the mode toward the outer edge and
        produces curvature (bend) radiation loss.  R=0 → straight (no-op).
        """
        if not self.bend_radius_um:
            return n2
        e = np.clip(2.0 * (x - self.bend_center_um) / self.bend_radius_um, -5.0, 5.0)
        return n2 * np.exp(e)[:, None]

    # -- index slice ----------------------------------------------------
    def index_slice(self, y: np.ndarray, bg_profile: np.ndarray,
                    bands: List[Tuple[Tuple[int, int], np.ndarray, float]],
                    layer_masks_at_z: Dict[Tuple[int, int], np.ndarray]) -> np.ndarray:
        """Build n(x, y) for one z slice.

        bg_profile : (ny,) background index vs y
        bands      : list of (mask_key, y-band bool (ny,), n_index)
        layer_masks_at_z : {mask_key: lateral bool mask (nx,) at this z}
        """
        nx = next(iter(layer_masks_at_z.values())).shape[0] if layer_masks_at_z else None
        if nx is None:
            return np.broadcast_to(bg_profile, (1, len(y))).copy()
        n_xy = np.broadcast_to(bg_profile[None, :], (nx, len(y))).copy()
        for key, band, n_val in bands:
            lat = layer_masks_at_z.get(key)
            if lat is None or not np.any(lat):
                continue
            # outer product of lateral mask (x) and vertical band (y)
            region = lat[:, None] & band[None, :]
            n_xy[region] = n_val
        return n_xy

    # -- input / output modes ------------------------------------------
    def gaussian_2d(self, x: np.ndarray, y: np.ndarray, x0: float, y0: float,
                    wx: float, wy: float) -> np.ndarray:
        gx = np.exp(-((x - x0) / max(wx, 1e-9)) ** 2)
        gy = np.exp(-((y - y0) / max(wy, 1e-9)) ** 2)
        field = (gx[:, None] * gy[None, :]).astype(np.complex128)
        return self.normalize(field, x, y)

    @staticmethod
    def normalize(field: np.ndarray, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        p = np.trapz(np.trapz(np.abs(field) ** 2, y, axis=1), x)
        if p <= 0:
            return field
        return field / np.sqrt(p)

    @staticmethod
    def power(field: np.ndarray, x: np.ndarray, y: np.ndarray) -> float:
        return float(np.trapz(np.trapz(np.abs(field) ** 2, y, axis=1), x))

    @staticmethod
    def overlap_power(field: np.ndarray, mode: np.ndarray,
                      x: np.ndarray, y: np.ndarray) -> float:
        ov = np.trapz(np.trapz(field * np.conj(mode), y, axis=1), x)
        return float(np.abs(ov) ** 2)

    # -- transverse operators (polarization-aware) ---------------------
    def _transverse_ops(self, x: np.ndarray, y: np.ndarray, alpha: complex,
                        n2: Optional[np.ndarray], pml_um: Optional[float] = None):
        """(opx, opy) tridiagonal transverse operators with prefactor ``alpha``
        and PML.  For TE the x-operator carries the index correction (E∥x);
        for TM the y-operator does (E∥y); scalar uses plain Laplacians.

        ``pml_um`` overrides the solver's PML width — mode solving passes 0 so
        the imaginary-distance iteration is not destabilised by the complex
        coordinate stretch (a guided mode is localised, so it needs no PML)."""
        p = self.pml_um if pml_um is None else pml_um
        s, m = self.pml_sigma, self.pml_m
        lx = lambda: _laplacian_tridiag(x, self.dx_um, alpha, p, s, m)
        ly = lambda: _laplacian_tridiag(y, self.dy_um, alpha, p, s, m)
        if self.pol == "TE":
            return (_corrected_laplacian_tridiag(x, self.dx_um, alpha, p, s, m, n2, axis=0), ly())
        if self.pol == "TM":
            return (lx(), _corrected_laplacian_tridiag(y, self.dy_um, alpha, p, s, m, n2, axis=1))
        return lx(), ly()

    def _transverse_ops_real(self, x: np.ndarray, y: np.ndarray, n2: np.ndarray):
        """Same operators with a real unit prefactor and NO PML — used only for
        the Rayleigh quotient that extracts n_eff during mode solving."""
        lx = lambda: _laplacian_tridiag(x, self.dx_um, 1.0, 0.0, 0.0, 1.0)
        ly = lambda: _laplacian_tridiag(y, self.dy_um, 1.0, 0.0, 0.0, 1.0)
        if self.pol == "TE":
            return (_corrected_laplacian_tridiag(x, self.dx_um, 1.0, 0.0, 0.0, 1.0, n2, axis=0), ly())
        if self.pol == "TM":
            return (lx(), _corrected_laplacian_tridiag(y, self.dy_um, 1.0, 0.0, 0.0, 1.0, n2, axis=1))
        return lx(), ly()

    @staticmethod
    def _adi_step(psi: np.ndarray, V: np.ndarray, opx, opy, dz: complex) -> np.ndarray:
        """One Peaceman–Rachford ADI step over (possibly complex) ``dz``.

        opx/opy = (lower, diag, upper) tridiag coefficients (1D or 2D); the
        index potential ``V`` is split symmetrically between the two half-steps.
        A real ``dz`` advances propagation; ``dz = -1j*dξ`` performs an
        imaginary-distance (mode-finding) step.
        """
        half = 0.5 * dz
        lx_low, lx_diag, lx_up = opx
        ly_low, ly_diag, ly_up = opy
        # half-step 1: implicit x, explicit y
        rhs = psi + half * (_tridiag_matvec(ly_low, ly_diag, ly_up, psi, axis=1) + 0.5 * V * psi)
        diag_x = 1.0 - half * (_as_axis_array(lx_diag, 0) + 0.5 * V)
        psi = _thomas_axis(-half * np.asarray(lx_low), diag_x, -half * np.asarray(lx_up), rhs, axis=0)
        # half-step 2: implicit y, explicit x
        rhs = psi + half * (_tridiag_matvec(lx_low, lx_diag, lx_up, psi, axis=0) + 0.5 * V * psi)
        diag_y = 1.0 - half * (_as_axis_array(ly_diag, 1) + 0.5 * V)
        psi = _thomas_axis(-half * np.asarray(ly_low), diag_y, -half * np.asarray(ly_up), rhs, axis=1)
        return psi

    def _rayleigh_neff(self, psi, n2, phys_ops, x, y) -> float:
        opx, opy = phys_ops
        Lpsi = (_tridiag_matvec(*opx, psi, axis=0) + _tridiag_matvec(*opy, psi, axis=1))
        num = np.trapz(np.trapz(np.conj(psi) * (Lpsi + self.k0 ** 2 * n2 * psi), y, axis=1), x)
        den = np.trapz(np.trapz(np.abs(psi) ** 2, y, axis=1), x)
        beta2 = float((num / den).real)
        return np.sqrt(max(beta2, 0.0)) / self.k0

    # -- eigenmode solver (imaginary-distance BPM) ---------------------
    def solve_mode(
        self,
        x: np.ndarray,
        y: np.ndarray,
        n2: np.ndarray,
        x0: float,
        y0: float,
        w0x: float,
        w0y: float,
        max_iter: int = 800,
        dxi_um: Optional[float] = None,
        tol: float = 1e-8,
        progress_cb: Optional[Callable[[float, str], None]] = None,
    ) -> Tuple[np.ndarray, float]:
        """Fundamental guided mode of the (z-invariant) cross-section ``n2``.

        Uses imaginary-distance ADI: the mode with the highest n_eff grows
        fastest, so iterate-and-renormalise converges to the fundamental.
        Returns (mode_field normalised, n_eff).  Polarization follows ``self.pol``.
        """
        alpha = 1j * (1.0 / (2.0 * self.beta_ref))
        if dxi_um is None:
            dxi_um = 4.0 * self.dz_um
        dz = -1j * float(dxi_um)
        n2 = self._bent(n2, x)        # bent-guide mode if a bend radius is set
        # No PML during mode solving: the complex coordinate stretch destabilises
        # the imaginary-distance iteration, and a guided mode is localised anyway.
        ops = self._transverse_ops(x, y, alpha, n2, pml_um=0.0)
        phys = self._transverse_ops_real(x, y, n2)
        V = alpha * (self.k0 ** 2 * n2 - self.beta_ref ** 2)

        psi = self.gaussian_2d(x, y, x0, y0, w0x, w0y)
        neff = self._rayleigh_neff(psi, n2, phys, x, y)
        neff_prev = None
        for it in range(max_iter):
            psi = self._adi_step(psi, V, ops[0], ops[1], dz)
            psi = self.normalize(psi, x, y)
            neff = self._rayleigh_neff(psi, n2, phys, x, y)
            if neff_prev is not None and abs(neff - neff_prev) < tol:
                break
            neff_prev = neff
            if progress_cb and it % 25 == 0:
                progress_cb(min(0.99, it / max_iter), f"Mode solve: iter {it}, n_eff={neff:.5f}")
        # rotate global phase so the field is (nearly) real-positive at its peak
        peak = psi.flat[int(np.argmax(np.abs(psi)))]
        if peak != 0:
            psi = psi * (np.conj(peak) / abs(peak))
        return self.normalize(psi, x, y), float(neff)

    # -- propagation ----------------------------------------------------
    def run(
        self,
        x: np.ndarray,
        y: np.ndarray,
        z: np.ndarray,
        index_fn: Callable[[int], np.ndarray],
        field0: np.ndarray,
        save_every: int = 4,
        progress_cb: Optional[Callable[[float, str], None]] = None,
        progress_prefix: str = "Running 3D BPM",
    ) -> BPM3DResult:
        """Propagate ``field0`` along z.

        index_fn(iz) -> n(x, y) for slice iz  (built lazily to save memory).
        """
        nx, ny, nz = len(x), len(y), len(z)
        dz = self.dz_um
        alpha = 1j * (1.0 / (2.0 * self.beta_ref))   # i * 1/(2 beta_ref)
        scalar = self.pol == "SCALAR"
        if scalar:
            opx, opy = self._transverse_ops(x, y, alpha, None)  # index-independent

        psi = field0.astype(np.complex128).copy()

        total_power = np.zeros(nz)
        total_power[0] = self.power(psi, x, y)
        topview = [np.trapz(np.abs(psi) ** 2, y, axis=1)]
        sideview = [np.trapz(np.abs(psi) ** 2, x, axis=0)]
        z_samples = [z[0]]

        stride = max(1, (nz - 1) // 200)
        save_every = max(1, int(save_every))
        if progress_cb:
            progress_cb(0.0, f"{progress_prefix} [{self.pol}]: nx={nx}, ny={ny}, nz={nz}")

        for iz in range(1, nz):
            n_xy = index_fn(iz - 1)
            n2 = self._bent(n_xy ** 2, x)
            V = alpha * (self.k0 ** 2 * n2 - self.beta_ref ** 2)
            if not scalar:
                opx, opy = self._transverse_ops(x, y, alpha, n2)   # rebuild per slice
            psi = self._adi_step(psi, V, opx, opy, dz)

            total_power[iz] = self.power(psi, x, y)
            if iz % save_every == 0 or iz == nz - 1:
                topview.append(np.trapz(np.abs(psi) ** 2, y, axis=1))
                sideview.append(np.trapz(np.abs(psi) ** 2, x, axis=0))
                z_samples.append(z[iz])

            if progress_cb and (iz % stride == 0 or iz == nz - 1):
                progress_cb((iz + 1) / nz, f"{progress_prefix}: step {iz + 1}/{nz}")

        return BPM3DResult(
            x=x, y=y, z_samples=np.array(z_samples),
            psi_out=psi, total_power=total_power,
            topview=np.array(topview), sideview=np.array(sideview),
            input_power=float(total_power[0]),
        )

    # -- full-vectorial propagation (Ex, Ey with cross coupling) --------
    def _cross_sources(self, Ex, Ey, ln_n2, alpha):
        """Off-diagonal coupling contributions to dE/dz (× the i/(2β) prefactor).

        Pxy Ey = d/dx[(∂ ln n² /∂y) Ey],  Pyx Ex = d/dy[(∂ ln n² /∂x) Ex].
        These are non-zero only where both transverse index gradients meet
        (i.e. at waveguide corners) — the physical origin of polarization mixing.
        """
        dln_dy = _ddy(ln_n2, self.dy_um)
        dln_dx = _ddx(ln_n2, self.dx_um)
        pxy = _ddx(dln_dy * Ey, self.dx_um)     # -> drives Ex
        pyx = _ddy(dln_dx * Ex, self.dy_um)     # -> drives Ey
        return alpha * pxy, alpha * pyx

    def _fullvec_diag_ops(self, x, y, n2, alpha, pml_um):
        """Diagonal transverse operators for full-vector: Ex is TE-like
        (correction along x), Ey is TM-like (correction along y)."""
        s, m = self.pml_sigma, self.pml_m
        opx_ex = _corrected_laplacian_tridiag(x, self.dx_um, alpha, pml_um, s, m, n2, axis=0)
        opy_ex = _laplacian_tridiag(y, self.dy_um, alpha, pml_um, s, m)
        opx_ey = _laplacian_tridiag(x, self.dx_um, alpha, pml_um, s, m)
        opy_ey = _corrected_laplacian_tridiag(y, self.dy_um, alpha, pml_um, s, m, n2, axis=1)
        return opx_ex, opy_ex, opx_ey, opy_ey

    def _fullvec_step(self, Ex, Ey, V, dops, ln_n2, alpha, dz):
        """One Strang-split full-vector step: cross ½ · diagonal ADI · cross ½."""
        opx_ex, opy_ex, opx_ey, opy_ey = dops
        sx, sy = self._cross_sources(Ex, Ey, ln_n2, alpha)
        Ex = Ex + 0.5 * dz * sx
        Ey = Ey + 0.5 * dz * sy
        Ex = self._adi_step(Ex, V, opx_ex, opy_ex, dz)
        Ey = self._adi_step(Ey, V, opx_ey, opy_ey, dz)
        sx, sy = self._cross_sources(Ex, Ey, ln_n2, alpha)
        Ex = Ex + 0.5 * dz * sx
        Ey = Ey + 0.5 * dz * sy
        return Ex, Ey

    def run_fullvec(
        self,
        x: np.ndarray,
        y: np.ndarray,
        z: np.ndarray,
        index_fn: Callable[[int], np.ndarray],
        ex0: np.ndarray,
        ey0: np.ndarray,
        save_every: int = 4,
        progress_cb: Optional[Callable[[float, str], None]] = None,
        progress_prefix: str = "Running full-vector 3D BPM",
    ) -> BPM3DResult:
        """Full-vectorial propagation of (Ex, Ey).

        Diagonal blocks Pxx/Pyy are the semi-vectorial TE (correction along x)
        and TM (correction along y) operators; the off-diagonal cross coupling
        is applied explicitly with Strang splitting (cross ½ · diagonal · cross ½),
        so polarization conversion at corners is captured while staying stable.
        """
        nx, ny, nz = len(x), len(y), len(z)
        dz = self.dz_um
        alpha = 1j * (1.0 / (2.0 * self.beta_ref))

        Ex = ex0.astype(np.complex128).copy()
        Ey = ey0.astype(np.complex128).copy()

        def tot_power(ex, ey):
            return (float(np.trapz(np.trapz(np.abs(ex) ** 2, y, axis=1), x))
                    + float(np.trapz(np.trapz(np.abs(ey) ** 2, y, axis=1), x)))

        total_power = np.zeros(nz)
        total_power[0] = tot_power(Ex, Ey)
        topv = lambda ex, ey: np.trapz(np.abs(ex) ** 2 + np.abs(ey) ** 2, y, axis=1)
        sidev = lambda ex, ey: np.trapz(np.abs(ex) ** 2 + np.abs(ey) ** 2, x, axis=0)
        topview = [topv(Ex, Ey)]
        sideview = [sidev(Ex, Ey)]
        z_samples = [z[0]]

        stride = max(1, (nz - 1) // 200)
        save_every = max(1, int(save_every))
        if progress_cb:
            progress_cb(0.0, f"{progress_prefix} [FULL]: nx={nx}, ny={ny}, nz={nz}")

        for iz in range(1, nz):
            n_xy = index_fn(iz - 1)
            n2 = self._bent(n_xy ** 2, x)
            ln_n2 = np.log(np.maximum(n2, 1e-12))
            V = alpha * (self.k0 ** 2 * n2 - self.beta_ref ** 2)
            dops = self._fullvec_diag_ops(x, y, n2, alpha, self.pml_um)
            Ex, Ey = self._fullvec_step(Ex, Ey, V, dops, ln_n2, alpha, dz)

            total_power[iz] = tot_power(Ex, Ey)
            if iz % save_every == 0 or iz == nz - 1:
                topview.append(topv(Ex, Ey))
                sideview.append(sidev(Ex, Ey))
                z_samples.append(z[iz])
            if progress_cb and (iz % stride == 0 or iz == nz - 1):
                progress_cb((iz + 1) / nz, f"{progress_prefix}: step {iz + 1}/{nz}")

        p_ex = float(np.trapz(np.trapz(np.abs(Ex) ** 2, y, axis=1), x))
        p_ey = float(np.trapz(np.trapz(np.abs(Ey) ** 2, y, axis=1), x))
        conv = p_ey / max(p_ex + p_ey, 1e-30)
        return BPM3DResult(
            x=x, y=y, z_samples=np.array(z_samples),
            psi_out=Ex, total_power=total_power,
            topview=np.array(topview), sideview=np.array(sideview),
            input_power=float(total_power[0]),
            psi_out_minor=Ey, conversion=conv,
        )

    # -- full-vector hybrid eigenmode solver ---------------------------
    @staticmethod
    def vector_overlap_power(ex, ey, mex, mey, x, y) -> float:
        """|⟨(Ex,Ey)|(mEx,mEy)⟩|² power overlap of two (normalised) vector fields."""
        ov = (np.trapz(np.trapz(ex * np.conj(mex), y, axis=1), x)
              + np.trapz(np.trapz(ey * np.conj(mey), y, axis=1), x))
        return float(np.abs(ov) ** 2)

    def _rayleigh_neff_fullvec(self, Ex, Ey, n2, ln_n2, real_ops, x, y) -> float:
        (oxe, oye, oxy, oyy) = real_ops
        dln_dy = _ddy(ln_n2, self.dy_um)
        dln_dx = _ddx(ln_n2, self.dx_um)
        HEx = (_tridiag_matvec(*oxe, Ex, axis=0) + _tridiag_matvec(*oye, Ex, axis=1)
               + _ddx(dln_dy * Ey, self.dx_um) + self.k0 ** 2 * n2 * Ex)
        HEy = (_tridiag_matvec(*oxy, Ey, axis=0) + _tridiag_matvec(*oyy, Ey, axis=1)
               + _ddy(dln_dx * Ex, self.dy_um) + self.k0 ** 2 * n2 * Ey)
        num = (np.trapz(np.trapz(np.conj(Ex) * HEx, y, axis=1), x)
               + np.trapz(np.trapz(np.conj(Ey) * HEy, y, axis=1), x))
        den = (np.trapz(np.trapz(np.abs(Ex) ** 2, y, axis=1), x)
               + np.trapz(np.trapz(np.abs(Ey) ** 2, y, axis=1), x))
        beta2 = float((num / den).real)
        return np.sqrt(max(beta2, 0.0)) / self.k0

    def solve_mode_fullvec(
        self, x, y, n2, x0, y0, w0x, w0y, launch="Ex",
        max_iter: int = 800, dxi_um: Optional[float] = None, tol: float = 1e-8,
        progress_cb: Optional[Callable[[float, str], None]] = None,
    ) -> Tuple[np.ndarray, np.ndarray, float]:
        """Fundamental **hybrid** mode (Ex, Ey) of cross-section ``n2`` via coupled
        imaginary-distance full-vector iteration.  Returns (Ex, Ey, n_eff).

        Seeds the launched component with a Gaussian; the cross coupling builds
        the minor component, so the converged field is the true vector mode."""
        alpha = 1j * (1.0 / (2.0 * self.beta_ref))
        if dxi_um is None:
            dxi_um = 4.0 * self.dz_um
        dz = -1j * float(dxi_um)
        n2 = self._bent(n2, x)        # bent-guide mode if a bend radius is set
        ln_n2 = np.log(np.maximum(n2, 1e-12))
        V = alpha * (self.k0 ** 2 * n2 - self.beta_ref ** 2)
        dops = self._fullvec_diag_ops(x, y, n2, alpha, pml_um=0.0)  # no PML (see solve_mode)
        # real, PML-free operators for the Rayleigh quotient
        real_ops = (
            _corrected_laplacian_tridiag(x, self.dx_um, 1.0, 0.0, 0.0, 1.0, n2, axis=0),
            _laplacian_tridiag(y, self.dy_um, 1.0, 0.0, 0.0, 1.0),
            _laplacian_tridiag(x, self.dx_um, 1.0, 0.0, 0.0, 1.0),
            _corrected_laplacian_tridiag(y, self.dy_um, 1.0, 0.0, 0.0, 1.0, n2, axis=1),
        )

        seed = self.gaussian_2d(x, y, x0, y0, w0x, w0y)
        Ex = seed.copy() if launch == "Ex" else np.zeros_like(seed)
        Ey = seed.copy() if launch == "Ey" else np.zeros_like(seed)

        def renorm(ex, ey):
            p = (np.trapz(np.trapz(np.abs(ex) ** 2, y, axis=1), x)
                 + np.trapz(np.trapz(np.abs(ey) ** 2, y, axis=1), x))
            if p > 0:
                s = 1.0 / np.sqrt(p)
                return ex * s, ey * s
            return ex, ey

        neff_prev = None
        neff = self._rayleigh_neff_fullvec(Ex, Ey, n2, ln_n2, real_ops, x, y)
        for it in range(max_iter):
            Ex, Ey = self._fullvec_step(Ex, Ey, V, dops, ln_n2, alpha, dz)
            Ex, Ey = renorm(Ex, Ey)
            neff = self._rayleigh_neff_fullvec(Ex, Ey, n2, ln_n2, real_ops, x, y)
            if neff_prev is not None and abs(neff - neff_prev) < tol:
                break
            neff_prev = neff
            if progress_cb and it % 25 == 0:
                progress_cb(min(0.99, it / max_iter), f"Vector mode: iter {it}, n_eff={neff:.5f}")
        # global-phase normalise on the dominant component's peak
        dom = Ex if launch == "Ex" else Ey
        peak = dom.flat[int(np.argmax(np.abs(dom)))]
        if peak != 0:
            ph = np.conj(peak) / abs(peak)
            Ex, Ey = Ex * ph, Ey * ph
        Ex, Ey = renorm(Ex, Ey)
        return Ex, Ey, float(neff)
