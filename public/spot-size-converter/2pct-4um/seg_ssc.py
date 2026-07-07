"""
Segmented (digitized) SSC inverse taper
=======================================

Reproduces and generalises the structure in ``silica_line_seg1.gds``:

  * two continuous lead-in blocks (a solid ridge, growing in width), then
  * a chain of *dashed* segments on a fixed pitch, where BOTH
      - the physical segment width  w(z)   grows  (tip -> bus), and
      - the duty cycle  d(z) = seg_len / pitch  shrinks (solid -> sparse)

The dashed segments act as a sub-wavelength/segmented waveguide: a low duty
lowers the local effective index, so the guide's confinement is ramped
gradually, converting the small waveguide mode to/from the large SMF mode.
Reverse direction (bus solid, tip dashed) is the usual "digital" SSC; this
file follows the uploaded convention (tip = wide-duty solid end on the left,
bus = narrow-duty sparse end... actually per the GDS the DUTY drops toward the
bus while WIDTH grows — both handled parametrically below).

Design parameters (defaults reproduce the uploaded GDS)
-------------------------------------------------------
  pitch          5.19 um        segment period (fixed)
  n_seg          60             number of dashed segments
  duty_start     0.79           duty at first dashed segment
  duty_end       0.25           duty at last dashed segment
  w_start        4.40 um        width at first dashed segment
  w_end          8.50 um        width at last dashed segment (= bus width)
  lead-ins       [(50,3.5),(40,4.3)]   (length_um, width_um) solid blocks

Interop: designed to plug into the same 3D-BPM evaluation as
``ssc_ga_optimizer.py`` (SCALAR/TE, 4.5 um silica core, 1301 nm).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import numpy as np


# ===========================================================================
# 1. Parametric segmented-SSC gene
# ===========================================================================

@dataclass
class SegSSCGene:
    """Digitized inverse taper defined by pitch + duty ramp + width ramp.

    All ramps are evaluated over the DASHED region only.  ``duty_profile`` and
    ``width_profile`` pick how the values interpolate between the start/end.
    """
    pitch_um: float = 5.19
    n_seg: int = 60
    duty_start: float = 0.79
    duty_end: float = 0.25
    w_start: float = 4.40
    w_end: float = 8.50
    duty_profile: str = "linear"      # "linear" | "quad" | "cos"
    width_profile: str = "linear"
    # continuous solid lead-ins at the tip end: list of (length_um, width_um)
    leadins: List[Tuple[float, float]] = field(
        default_factory=lambda: [(50.0, 3.5), (40.0, 4.3)])
    y_center_um: float = 0.0

    # ---- ramp helpers ----
    @staticmethod
    def _ramp(a: float, b: float, t: np.ndarray, kind: str) -> np.ndarray:
        if kind == "quad":
            t = t ** 2
        elif kind == "cos":
            t = 0.5 * (1 - np.cos(np.pi * t))
        return a + (b - a) * t

    def seg_index_t(self) -> np.ndarray:
        return np.linspace(0.0, 1.0, self.n_seg)

    def duties(self) -> np.ndarray:
        return self._ramp(self.duty_start, self.duty_end,
                          self.seg_index_t(), self.duty_profile)

    def widths(self) -> np.ndarray:
        return self._ramp(self.w_start, self.w_end,
                         self.seg_index_t(), self.width_profile)

    @property
    def leadin_length(self) -> float:
        return float(sum(L for L, _ in self.leadins))

    @property
    def dashed_length(self) -> float:
        return self.n_seg * self.pitch_um

    @property
    def total_length(self) -> float:
        return self.leadin_length + self.dashed_length

    # ---- segment rectangles (in local z starting at 0) ----
    def segments(self) -> List[Tuple[float, float, float]]:
        """Return list of (z_start, z_end, width) for every SOLID block:
        the lead-ins first, then each dashed tooth (duty-limited length)."""
        rects = []
        z = 0.0
        for L, w in self.leadins:
            rects.append((z, z + L, w))
            z += L
        duties = self.duties()
        widths = self.widths()
        for i in range(self.n_seg):
            seg_len = duties[i] * self.pitch_um
            z0 = z + 0.5 * (self.pitch_um - seg_len)   # centre the tooth in its period
            rects.append((z0, z0 + seg_len, widths[i]))
            z += self.pitch_um
        return rects


# ===========================================================================
# 2. GDS export (matches the uploaded style)
# ===========================================================================

def seg_gene_to_gds(gene: SegSSCGene, path: str, layer: int = 0, datatype: int = 0,
                    cell_name: str = "CIRCUIT", y_offset_um: float = 0.0) -> str:
    """Write the segmented SSC to GDS. One rectangle per solid block."""
    import gdstk
    lib = gdstk.Library()
    cell = lib.new_cell(cell_name)
    yc = gene.y_center_um + y_offset_um
    for (z0, z1, w) in gene.segments():
        hw = 0.5 * w
        cell.add(gdstk.rectangle((z0, yc - hw), (z1, yc + hw),
                                 layer=layer, datatype=datatype))
    lib.write_gds(path)
    return path


# ===========================================================================
# 3. BPM index field for the segmented SSC
# ===========================================================================

class SegSSCIndexField:
    """Lazy n(x, y; z) for the digitized inverse taper.

    Vertical: buried silica core of fixed thickness (tanh soft edges).
    Lateral+longitudinal: a z-slice is CORE only where it falls inside a solid
    block (between teeth -> pure cladding).  A short longitudinal soft edge
    (~pitch/8) rounds the tooth ends so the BPM sees no hard index step.
    """

    def __init__(self, gene: SegSSCGene, n_core: float, n_clad: float,
                 core_thick_um: float, x: np.ndarray, y: np.ndarray,
                 z: np.ndarray, edge_um: float = 0.12,
                 z_edge_um: Optional[float] = None):
        self.g = gene
        self.x = x
        self.y = y
        self.z = z
        self.edge = edge_um
        self.z_edge = z_edge_um if z_edge_um is not None else gene.pitch_um / 8.0
        self._nclad = n_clad
        self._dn = n_core - n_clad
        ht = 0.5 * core_thick_um
        yc = gene.y_center_um
        self._vert = 0.5 * (np.tanh((y - yc + ht) / edge_um)
                            - np.tanh((y - yc - ht) / edge_um))     # (ny,)
        self._rects = gene.segments()
        # precompute, for each z sample, the active (width, longitudinal weight)
        self._zw = self._build_zslice_widths(z)

    def _build_zslice_widths(self, z: np.ndarray) -> np.ndarray:
        """For each z, effective (width, on_weight) from the nearest tooth.

        Vectorised: build the smooth top-hat membership of every rect over the
        whole z grid at once, take the per-z max, and pick the width of the
        argmax rect."""
        rects = np.array(self._rects)          # (Nrect, 3): z0, z1, w
        z0 = rects[:, 0][:, None]              # (Nrect, 1)
        z1 = rects[:, 1][:, None]
        ww = rects[:, 2]                       # (Nrect,)
        zz = z[None, :]                        # (1, nz)
        m = 0.5 * (np.tanh((zz - z0) / self.z_edge)
                   - np.tanh((zz - z1) / self.z_edge))       # (Nrect, nz)
        best = np.argmax(m, axis=0)            # (nz,) index of dominant rect
        on = m[best, np.arange(len(z))]        # (nz,) its membership weight
        wsel = ww[best]                        # (nz,) its width
        return np.column_stack([wsel, on])

    def __call__(self, iz: int) -> np.ndarray:
        w, on = self._zw[iz]
        hw = 0.5 * w
        lat = 0.5 * (np.tanh((self.x + hw) / self.edge)
                     - np.tanh((self.x - hw) / self.edge))          # (nx,)
        mask = on * np.outer(lat, self._vert)                       # (nx, ny)
        return self._nclad + self._dn * mask


# ===========================================================================
# 4. Convenience: build the default (uploaded-matching) gene
# ===========================================================================

def default_uploaded_gene() -> SegSSCGene:
    return SegSSCGene(
        pitch_um=5.19, n_seg=60,
        duty_start=0.79, duty_end=0.25,
        w_start=4.40, w_end=8.50,
        leadins=[(50.0, 3.5), (40.0, 4.3)],
    )


# ===========================================================================
# 5. 3D-BPM evaluation  (SMF at the sparse/wide "bus" facet)
# ===========================================================================

def evaluate_seg_gene(gene: SegSSCGene, plat, grid,
                      progress_cb=None) -> dict:
    """Coupling loss [dB] of a segmented SSC at plat.lam_um.

    Facet convention (matches the uploaded GDS): the LOW-duty, wide-width end
    (last dashed segment) is the low-index SMF facet; the solid, narrow lead-in
    end is the chip waveguide.  We launch the SMF-28 Gaussian at the SMF facet
    and overlap the chip guided mode at the opposite facet (reciprocity gives
    the same number as chip->SMF)."""
    from ssc_ga_optimizer import BPM3DSolver

    dx, dy = grid.dx_um, grid.dy_um
    x = np.arange(-grid.x_half_um, grid.x_half_um + 1e-9, dx)
    y = np.arange(-grid.y_half_um, grid.y_half_um + 1e-9, dy)
    L = gene.total_length
    nz = max(4, int(round(L / grid.dz_um)) + 1)
    z = np.linspace(0.0, L, nz)

    idx = SegSSCIndexField(gene, plat.n_core, plat.n_clad,
                           plat.core_thick_um, x, y, z)
    n_ref = 0.5 * (plat.n_core + plat.n_clad)
    solver = BPM3DSolver(wavelength_um=plat.lam_um, n_ref=n_ref,
                         dx_um=dx, dy_um=dy, dz_um=(z[1] - z[0]),
                         pml_um=grid.pml_um, pol=grid.pol)

    # launch at the SMF facet (last dashed segment = idx(nz-1)); propagate to chip
    rev = lambda iz: idx(nz - 1 - iz)
    w0 = plat.smf_w0_um
    f0 = solver.normalize(solver.gaussian_2d(x, y, 0.0, 0.0, w0, w0), x, y)
    res = solver.run(x, y, z, rev, f0, save_every=grid.save_every,
                     progress_cb=progress_cb, progress_prefix="SegSSC BPM")

    n_chip = idx(0)                    # solid narrow chip cross-section
    chip_mode, neff = solver.solve_mode(
        x, y, n_chip ** 2, 0.0, 0.0, 2.0, plat.core_thick_um * 0.5)
    eta = float(np.clip(solver.overlap_power(res.psi_out, chip_mode, x, y),
                        1e-12, 1.0))
    return dict(coupling_loss_dB=-10.0 * np.log10(eta), efficiency=eta,
                chip_neff=neff,
                residual_power=float(res.total_power[-1] / max(res.total_power[0], 1e-30)),
                length_um=L, result=res, x=x, y=y, z=z)


# ===========================================================================
# 6. Genetic algorithm over the segmented-SSC parameters
# ===========================================================================

from dataclasses import dataclass as _dc


@_dc
class SegGABounds:
    pitch_um: Tuple[float, float] = (3.0, 7.0)
    duty_start: Tuple[float, float] = (0.6, 0.95)
    duty_end: Tuple[float, float] = (0.12, 0.45)
    w_start: Tuple[float, float] = (3.0, 5.5)
    w_end: float = 8.5                       # fixed bus width
    n_seg: Tuple[int, int] = (40, 90)
    profiles: Tuple[str, ...] = ("linear", "quad", "cos")


class SegSSCGeneticOptimizer:
    def __init__(self, plat, grid, bounds: SegGABounds,
                 pop_size=16, n_gen=12, elite=3, mut_rate=0.3,
                 seed=0, leadins=None, log_cb=None):
        self.plat = plat
        self.grid = grid
        self.b = bounds
        self.pop_size = pop_size
        self.n_gen = n_gen
        self.elite = elite
        self.mut_rate = mut_rate
        self.rng = np.random.default_rng(seed)
        self.leadins = leadins if leadins is not None else [(50.0, 3.5), (40.0, 4.3)]
        self.log = log_cb or (lambda s: print(s, flush=True))
        self.history = []

    def _rand(self) -> SegSSCGene:
        b = self.b
        return SegSSCGene(
            pitch_um=self.rng.uniform(*b.pitch_um),
            n_seg=int(self.rng.integers(b.n_seg[0], b.n_seg[1] + 1)),
            duty_start=self.rng.uniform(*b.duty_start),
            duty_end=self.rng.uniform(*b.duty_end),
            w_start=self.rng.uniform(*b.w_start),
            w_end=b.w_end,
            duty_profile=self.rng.choice(b.profiles),
            width_profile=self.rng.choice(b.profiles),
            leadins=list(self.leadins),
        )

    def _clip(self, g: SegSSCGene) -> SegSSCGene:
        b = self.b
        g.pitch_um = float(np.clip(g.pitch_um, *b.pitch_um))
        g.n_seg = int(np.clip(g.n_seg, *b.n_seg))
        g.duty_start = float(np.clip(g.duty_start, *b.duty_start))
        g.duty_end = float(np.clip(g.duty_end, *b.duty_end))
        g.w_start = float(np.clip(g.w_start, *b.w_start))
        g.w_end = b.w_end
        return g

    def _cross(self, a, b_):
        w = self.rng.random()
        mix = lambda pa, pb: w * pa + (1 - w) * pb
        return self._clip(SegSSCGene(
            pitch_um=mix(a.pitch_um, b_.pitch_um),
            n_seg=int(round(mix(a.n_seg, b_.n_seg))),
            duty_start=mix(a.duty_start, b_.duty_start),
            duty_end=mix(a.duty_end, b_.duty_end),
            w_start=mix(a.w_start, b_.w_start),
            w_end=self.b.w_end,
            duty_profile=a.duty_profile if self.rng.random() < 0.5 else b_.duty_profile,
            width_profile=a.width_profile if self.rng.random() < 0.5 else b_.width_profile,
            leadins=list(self.leadins),
        ))

    def _mutate(self, g):
        r = self.rng
        if r.random() < self.mut_rate: g.pitch_um += r.normal(0, 0.5)
        if r.random() < self.mut_rate: g.n_seg += int(r.normal(0, 6))
        if r.random() < self.mut_rate: g.duty_start += r.normal(0, 0.06)
        if r.random() < self.mut_rate: g.duty_end += r.normal(0, 0.05)
        if r.random() < self.mut_rate: g.w_start += r.normal(0, 0.3)
        if r.random() < self.mut_rate: g.duty_profile = r.choice(self.b.profiles)
        if r.random() < self.mut_rate: g.width_profile = r.choice(self.b.profiles)
        return self._clip(g)

    def _fit(self, g):
        try:
            return evaluate_seg_gene(g, self.plat, self.grid)
        except Exception as e:
            return dict(coupling_loss_dB=99.0, efficiency=0.0, error=str(e),
                        length_um=g.total_length)

    def run(self):
        import time
        pop = [self._rand() for _ in range(self.pop_size)]
        best = None
        t0 = time.time()
        for gen in range(self.n_gen):
            scored = []
            for g in pop:
                f = self._fit(g)
                scored.append((f["coupling_loss_dB"], g, f))
            scored.sort(key=lambda t: t[0])
            gb_loss, gb_gene, gb = scored[0]
            if best is None or gb_loss < best[0]:
                best = (gb_loss, gb_gene, gb)
            self.history.append(dict(gen=gen, best_loss=best[0]))
            self.log(f"[gen {gen:02d}] best={best[0]:.4f} dB  "
                     f"(eff={best[2]['efficiency']*100:.2f}%, "
                     f"L={best[1].total_length:.0f} um, "
                     f"pitch={best[1].pitch_um:.2f}, nseg={best[1].n_seg})  "
                     f"[{time.time()-t0:.0f}s]")
            elites = [s[1] for s in scored[:self.elite]]
            parents = [s[1] for s in scored[:max(self.elite, self.pop_size // 2)]]
            newpop = list(elites)
            while len(newpop) < self.pop_size:
                a, c = self.rng.choice(len(parents), 2, replace=True)
                newpop.append(self._mutate(self._cross(parents[a], parents[c])))
            pop = newpop
        return dict(best_loss_dB=best[0], best_gene=best[1], best_eval=best[2],
                    history=self.history)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Segmented (digitized) SSC generator / optimizer")
    ap.add_argument("--gds", default="seg_ssc.gds")
    ap.add_argument("--pitch", type=float, default=5.19)
    ap.add_argument("--nseg", type=int, default=60)
    ap.add_argument("--duty0", type=float, default=0.79)
    ap.add_argument("--duty1", type=float, default=0.25)
    ap.add_argument("--w0", type=float, default=4.40)
    ap.add_argument("--w1", type=float, default=8.50)
    ap.add_argument("--optimize", action="store_true", help="run GA optimization")
    ap.add_argument("--evaluate", action="store_true", help="3D-BPM evaluate the current gene")
    ap.add_argument("--pop", type=int, default=14)
    ap.add_argument("--gen", type=int, default=10)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--fast", action="store_true")
    args = ap.parse_args()

    if args.optimize or args.evaluate:
        from ssc_ga_optimizer import Platform, GridConfig
        plat = Platform()
        grid = GridConfig(dx_um=0.18, dy_um=0.18, dz_um=0.4,
                          x_half_um=18.0, y_half_um=15.0, pml_um=3.5, save_every=15)
        if args.fast:
            grid = GridConfig(dx_um=0.25, dy_um=0.25, dz_um=0.6,
                              x_half_um=16.0, y_half_um=13.0, pml_um=3.0, save_every=20)
        print(f"Platform: n_clad={plat.n_clad:.5f} n_core={plat.n_core:.5f} "
              f"t={plat.core_thick_um}um SMF w0={plat.smf_w0_um:.2f}um")

    if args.optimize:
        ga = SegSSCGeneticOptimizer(plat, grid, SegGABounds(),
                                    pop_size=args.pop, n_gen=args.gen, seed=args.seed)
        result = ga.run()
        g = result["best_gene"]
        print("\n=== BEST SEGMENTED SSC ===")
        print(f"coupling loss @1301nm : {result['best_loss_dB']:.4f} dB")
        print(f"efficiency            : {result['best_eval']['efficiency']*100:.3f} %")
        print(f"pitch / n_seg         : {g.pitch_um:.3f} um / {g.n_seg}")
        print(f"duty  start->end      : {g.duty_start:.3f} -> {g.duty_end:.3f} ({g.duty_profile})")
        print(f"width start->end      : {g.w_start:.3f} -> {g.w_end:.3f} um ({g.width_profile})")
        print(f"total length          : {g.total_length:.1f} um")
        seg_gene_to_gds(g, args.gds, y_offset_um=200.0)
        print(f"GDS written           : {args.gds}")
    elif args.evaluate:
        g = SegSSCGene(pitch_um=args.pitch, n_seg=args.nseg,
                       duty_start=args.duty0, duty_end=args.duty1,
                       w_start=args.w0, w_end=args.w1)
        ev = evaluate_seg_gene(g, plat, grid)
        print(f"coupling loss @1301nm : {ev['coupling_loss_dB']:.4f} dB")
        print(f"efficiency            : {ev['efficiency']*100:.3f} %")
        print(f"chip n_eff / residual : {ev['chip_neff']:.5f} / {ev['residual_power']:.3f}")
        seg_gene_to_gds(g, args.gds, y_offset_um=200.0)
        print(f"GDS written           : {args.gds}")
    else:
        g = SegSSCGene(pitch_um=args.pitch, n_seg=args.nseg,
                       duty_start=args.duty0, duty_end=args.duty1,
                       w_start=args.w0, w_end=args.w1)
        p = seg_gene_to_gds(g, args.gds, y_offset_um=200.0)
        print(f"segments (incl. lead-ins): {len(g.segments())}")
        print(f"dashed length : {g.dashed_length:.1f} um")
        print(f"total length  : {g.total_length:.1f} um")
        print(f"GDS written   : {p}")
