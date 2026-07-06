"""
SSC (Spot-Size Converter) inverse-taper GA optimizer
=====================================================

Goal
----
Find the segmented lateral inverse-taper that minimises the SMF-28 <-> chip
coupling loss at 1301 nm, on a 2% Delta silica platform with a 4.5 um thick
core (thickness fixed; only the lateral width w(z) varies).

Pipeline per GA candidate
--------------------------
  gene  ->  width profile w(z)  ->  (a) GDS polygon (saved for the best)
                                    (b) index_fn n(x,y;z) for the 3D BPM
  launch SMF-28 Gaussian at the WIDE (chip) end going toward the TIP,
  overlap the output at the tip with the SMF-28 mode  ->  coupling efficiency.

Reciprocity note
----------------
Coupling SMF->chip == chip->SMF (linear, lossless-medium reciprocity), so we
launch the SMF Gaussian and overlap against the SMF Gaussian at the opposite
facet.  The taper does the mode-size conversion in between.  Loss reported is
the mode-overlap (mismatch) loss of the taper; propagation/radiation loss is
captured by the BPM power that leaks out of the overlap.

Platform (computed at 1301 nm, fused-silica Sellmeier)
------------------------------------------------------
  n_clad = 1.44691 ,  n_core(2%) = 1.47674 ,  core thickness = 4.5 um
  SMF-28 MFD @1310nm ~ 9.2 um  ->  Gaussian waist w0 = 4.6 um

Dependencies: numpy, gdstk, and the local bpm3d.py (BPM3DSolver).
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Callable, List, Optional, Tuple

import numpy as np

from bpm3d import BPM3DSolver


# ===========================================================================
# 1. Platform constants
# ===========================================================================

def fused_silica_index(lam_um: float) -> float:
    """Sellmeier refractive index of fused silica."""
    B = (0.6961663, 0.4079426, 0.8974794)
    C = (0.0684043 ** 2, 0.1162414 ** 2, 9.896161 ** 2)
    l2 = lam_um ** 2
    n2 = 1.0 + sum(B[i] * l2 / (l2 - C[i]) for i in range(3))
    return float(np.sqrt(n2))


@dataclass
class Platform:
    lam_um: float = 1.301
    delta: float = 0.02                # 2 %
    core_thick_um: float = 4.5         # fixed vertical thickness
    smf_mfd_um: float = 9.2            # SMF-28 mode-field diameter @1310nm

    n_clad: float = field(init=False)
    n_core: float = field(init=False)
    smf_w0_um: float = field(init=False)

    def __post_init__(self):
        self.n_clad = fused_silica_index(self.lam_um)
        self.n_core = self.n_clad / np.sqrt(1.0 - 2.0 * self.delta)
        self.smf_w0_um = 0.5 * self.smf_mfd_um


# ===========================================================================
# 2. Taper geometry (gene -> width profile)
# ===========================================================================

@dataclass
class TaperGene:
    """A segmented lateral inverse taper described by breakpoints.

    tip_w   : width at the SMF (input) facet          [um]
    seg_L   : length of each segment                  [um]  (len = n_seg)
    seg_w   : width at the END of each segment        [um]  (len = n_seg)

    Node widths are  [tip_w, seg_w[0], seg_w[1], ... seg_w[-1]]  at the
    cumulative-length nodes  [0, L0, L0+L1, ...].  The last node width is the
    bus-waveguide width (the wide, chip end).
    """
    tip_w: float
    seg_L: List[float]
    seg_w: List[float]

    @property
    def total_length(self) -> float:
        return float(sum(self.seg_L))

    def node_positions(self) -> np.ndarray:
        return np.concatenate([[0.0], np.cumsum(self.seg_L)])

    def node_widths(self) -> np.ndarray:
        return np.concatenate([[self.tip_w], np.asarray(self.seg_w, float)])

    def width_at(self, z: np.ndarray) -> np.ndarray:
        """Piecewise-linear half-width interpolation w(z) over [0, L]."""
        zp = self.node_positions()
        wp = self.node_widths()
        return np.interp(np.clip(z, 0.0, zp[-1]), zp, wp)


# ===========================================================================
# 3. GDS export (for the winning geometry / mask)
# ===========================================================================

def gene_to_gds(gene: TaperGene, path: str, layer: int = 1, datatype: int = 0,
                n_pts: int = 400, bus_extra_um: float = 20.0,
                tip_extra_um: float = 10.0, cell_name: str = "SSC_TAPER") -> str:
    """Write the taper (plus straight lead-ins) to a GDS file. Returns path."""
    import gdstk

    L = gene.total_length
    z = np.linspace(0.0, L, n_pts)
    hw = 0.5 * gene.width_at(z)                       # half-width

    # add straight tip lead-in (z<0) and bus lead-out (z>L)
    z_full = np.concatenate([[-tip_extra_um], z, [L + bus_extra_um]])
    hw_full = np.concatenate([[hw[0]], hw, [hw[-1]]])

    top = np.column_stack([z_full, hw_full])
    bot = np.column_stack([z_full[::-1], -hw_full[::-1]])
    pts = np.vstack([top, bot])

    lib = gdstk.Library()
    cell = lib.new_cell(cell_name)
    cell.add(gdstk.Polygon(pts, layer=layer, datatype=datatype))
    lib.write_gds(path)
    return path


# ===========================================================================
# 4. Index builder for the 3D BPM  (thickness fixed = buried silica core)
# ===========================================================================

class TaperIndexField:
    """Lazy n(x, y; z) for a lateral inverse taper with fixed vertical core.

    Vertical: a single buried core slab of thickness ``core_thick_um`` centred
    on y = 0, index n_core; elsewhere n_clad (symmetric buried silica).
    Lateral: |x| < half-width(z) -> core.  A tanh soft edge (~0.15 um) removes
    the index-staircase noise so the BPM sees a smooth guide.
    """

    def __init__(self, gene: TaperGene, plat: Platform,
                 x: np.ndarray, y: np.ndarray, z: np.ndarray,
                 edge_um: float = 0.12):
        self.gene = gene
        self.plat = plat
        self.x = x
        self.y = y
        self.z = z
        self.edge = edge_um
        # vertical core mask (z-invariant)  -> shape (ny,)
        ht = 0.5 * plat.core_thick_um
        self._vert = 0.5 * (np.tanh((y + ht) / edge_um) - np.tanh((y - ht) / edge_um))
        self._hw = 0.5 * gene.width_at(z)             # half-widths per slice
        self._dn = plat.n_core - plat.n_clad
        self._nclad = plat.n_clad

    def __call__(self, iz: int) -> np.ndarray:
        hw = self._hw[iz]
        lat = 0.5 * (np.tanh((self.x + hw) / self.edge)
                     - np.tanh((self.x - hw) / self.edge))   # (nx,)
        mask = np.outer(lat, self._vert)                     # (nx, ny) in [0,1]
        return self._nclad + self._dn * mask


# ===========================================================================
# 5. Single-candidate evaluation  ->  coupling loss [dB]
# ===========================================================================

@dataclass
class GridConfig:
    dx_um: float = 0.15
    dy_um: float = 0.15
    dz_um: float = 0.25
    x_half_um: float = 14.0
    y_half_um: float = 12.0
    pml_um: float = 2.5
    pol: str = "SCALAR"          # SCALAR is plenty for a symmetric silica taper
    save_every: int = 8


def _make_grids(gene: TaperGene, g: GridConfig):
    x = np.arange(-g.x_half_um, g.x_half_um + 1e-9, g.dx_um)
    y = np.arange(-g.y_half_um, g.y_half_um + 1e-9, g.dy_um)
    L = gene.total_length
    nz = max(4, int(round(L / g.dz_um)) + 1)
    z = np.linspace(0.0, L, nz)
    return x, y, z


def evaluate_gene(gene: TaperGene, plat: Platform, g: GridConfig,
                  progress_cb: Optional[Callable] = None) -> dict:
    """Return dict with coupling_loss_dB, efficiency, and diagnostics."""
    x, y, z = _make_grids(gene, g)

    n_ref = 0.5 * (plat.n_core + plat.n_clad)
    solver = BPM3DSolver(
        wavelength_um=plat.lam_um, n_ref=n_ref,
        dx_um=g.dx_um, dy_um=g.dy_um, dz_um=(z[1] - z[0]),
        pml_um=g.pml_um, pol=g.pol,
    )

    idx = TaperIndexField(gene, plat, x, y, z)

    # Inverse taper: the NARROW tip (z=0) faces the SMF; the WIDE bus (z=L)
    # is the chip waveguide.  Launch the SMF-28 Gaussian at the tip and
    # propagate toward the bus; at the bus facet overlap against the guided
    # chip mode of the bus cross-section (solved once).  Reciprocity makes
    # this equal to the chip->SMF direction.
    nz = len(z)

    w0 = plat.smf_w0_um
    field0 = solver.gaussian_2d(x, y, 0.0, 0.0, w0, w0)     # circular SMF mode
    field0 = solver.normalize(field0, x, y)

    res = solver.run(x, y, z, idx, field0, save_every=g.save_every,
                     progress_cb=progress_cb, progress_prefix="SSC BPM")

    # Guided mode of the wide bus cross-section (chip waveguide) at z = L.
    n_bus = idx(nz - 1)
    bus_hw = 0.5 * gene.node_widths()[-1]
    seed_w = max(bus_hw, plat.core_thick_um * 0.5)
    bus_mode, _neff = solver.solve_mode(
        x, y, n_bus ** 2, 0.0, 0.0, seed_w, plat.core_thick_um * 0.5,
    )
    eta_overlap = solver.overlap_power(res.psi_out, bus_mode, x, y)   # mode match
    eff = float(np.clip(eta_overlap, 1e-12, 1.0))
    loss_dB = -10.0 * np.log10(eff)

    return dict(
        coupling_loss_dB=loss_dB,
        efficiency=eff,
        residual_power=float(res.total_power[-1] / max(res.total_power[0], 1e-30)),
        length_um=gene.total_length,
        result=res, x=x, y=y, z=z,
    )


# ===========================================================================
# 6. Genetic algorithm
# ===========================================================================

@dataclass
class GABounds:
    tip_w: Tuple[float, float] = (0.20, 1.20)     # inverse-taper tip width [um]
    seg_w: Tuple[float, float] = (0.5, 8.0)       # intermediate/end widths [um]
    seg_L: Tuple[float, float] = (30.0, 400.0)    # per-segment length [um]
    bus_w: float = 4.5                            # fixed wide-end (chip) width
    n_seg: int = 3                                # number of taper segments


class SSCGeneticOptimizer:
    def __init__(self, plat: Platform, grid: GridConfig, bounds: GABounds,
                 pop_size: int = 16, n_gen: int = 12,
                 elite: int = 3, mut_rate: float = 0.25,
                 seed: Optional[int] = 0,
                 log_cb: Optional[Callable[[str], None]] = None):
        self.plat = plat
        self.grid = grid
        self.b = bounds
        self.pop_size = pop_size
        self.n_gen = n_gen
        self.elite = elite
        self.mut_rate = mut_rate
        self.rng = np.random.default_rng(seed)
        self.log = log_cb or (lambda s: print(s, flush=True))
        self.history: List[dict] = []

    # ---- gene <-> vector ----
    def _random_gene(self) -> TaperGene:
        b = self.b
        tip = self.rng.uniform(*b.tip_w)
        Ls = self.rng.uniform(b.seg_L[0], b.seg_L[1], b.n_seg).tolist()
        # monotonic increasing widths from tip -> bus (enforce adiabatic-ish)
        inner = np.sort(self.rng.uniform(b.seg_w[0], b.bus_w, b.n_seg - 1))
        ws = inner.tolist() + [b.bus_w]
        return TaperGene(tip_w=tip, seg_L=Ls, seg_w=ws)

    def _clip_gene(self, g: TaperGene) -> TaperGene:
        b = self.b
        tip = float(np.clip(g.tip_w, *b.tip_w))
        Ls = [float(np.clip(L, *b.seg_L)) for L in g.seg_L]
        ws = [float(np.clip(w, b.seg_w[0], b.bus_w)) for w in g.seg_w[:-1]]
        ws = sorted(ws) + [b.bus_w]                # keep monotonic + fixed bus
        # enforce tip <= first inner width
        if ws and tip > ws[0]:
            tip = min(tip, ws[0])
        return TaperGene(tip_w=tip, seg_L=Ls, seg_w=ws)

    def _crossover(self, a: TaperGene, b: TaperGene) -> TaperGene:
        w = self.rng.random()
        tip = w * a.tip_w + (1 - w) * b.tip_w
        Ls = [w * la + (1 - w) * lb for la, lb in zip(a.seg_L, b.seg_L)]
        ws = [w * wa + (1 - w) * wb for wa, wb in zip(a.seg_w, b.seg_w)]
        return self._clip_gene(TaperGene(tip, Ls, ws))

    def _mutate(self, g: TaperGene) -> TaperGene:
        b = self.b
        tip = g.tip_w
        Ls = list(g.seg_L)
        ws = list(g.seg_w)
        if self.rng.random() < self.mut_rate:
            tip += self.rng.normal(0, 0.15)
        for i in range(len(Ls)):
            if self.rng.random() < self.mut_rate:
                Ls[i] += self.rng.normal(0, 40.0)
        for i in range(len(ws) - 1):              # never mutate the fixed bus width
            if self.rng.random() < self.mut_rate:
                ws[i] += self.rng.normal(0, 0.6)
        return self._clip_gene(TaperGene(tip, Ls, ws))

    # ---- fitness ----
    def _fitness(self, g: TaperGene) -> dict:
        try:
            out = evaluate_gene(g, self.plat, self.grid)
        except Exception as e:                    # a broken gene shouldn't kill the run
            return dict(coupling_loss_dB=99.0, efficiency=0.0, error=str(e),
                        length_um=g.total_length)
        return out

    def run(self) -> dict:
        pop = [self._random_gene() for _ in range(self.pop_size)]
        best = None
        t0 = time.time()
        for gen in range(self.n_gen):
            scored = []
            for i, g in enumerate(pop):
                f = self._fitness(g)
                scored.append((f["coupling_loss_dB"], g, f))
            scored.sort(key=lambda t: t[0])
            gen_best = scored[0]
            if best is None or gen_best[0] < best[0]:
                best = gen_best
            self.history.append(dict(gen=gen, best_loss=gen_best[0],
                                     best_len=gen_best[1].total_length))
            self.log(f"[gen {gen:02d}] best={gen_best[0]:.4f} dB  "
                     f"(eff={gen_best[2]['efficiency']*100:.2f}%, "
                     f"L={gen_best[1].total_length:.0f} um)  "
                     f"[{time.time()-t0:.0f}s]")

            # next generation: elitism + crossover/mutation of tournament winners
            elites = [s[1] for s in scored[:self.elite]]
            newpop = list(elites)
            parents = [s[1] for s in scored[:max(self.elite, self.pop_size // 2)]]
            while len(newpop) < self.pop_size:
                a, b = self.rng.choice(len(parents), 2, replace=True)
                child = self._crossover(parents[a], parents[b])
                child = self._mutate(child)
                newpop.append(child)
            pop = newpop

        return dict(best_loss_dB=best[0], best_gene=best[1], best_eval=best[2],
                    history=self.history)


# ===========================================================================
# 7. Result visualization
# ===========================================================================

def plot_best(gene: TaperGene, ev: dict, plat: Platform, path: str) -> str:
    """Save a 4-panel summary: width profile, side/top BPM views, output field."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    res = ev["result"]
    x, y, z = ev["x"], ev["y"], ev["z"]
    fig, ax = plt.subplots(2, 2, figsize=(11, 7))

    # (a) width profile
    zz = np.linspace(0, gene.total_length, 400)
    ax[0, 0].plot(zz, gene.width_at(zz), "b-")
    ax[0, 0].plot(gene.node_positions(), gene.node_widths(), "ro", ms=4)
    ax[0, 0].set_xlabel("z (propagation) [um]")
    ax[0, 0].set_ylabel("Core width [um]")
    ax[0, 0].set_title(f"Inverse-taper profile  (tip={gene.tip_w:.2f} um, L={gene.total_length:.0f} um)")
    ax[0, 0].grid(alpha=0.3)

    # (b) top view |E|^2 integrated over y  (z vs x)
    tv = res.topview.T                       # (nx, nsaved)
    ax[0, 1].imshow(tv, aspect="auto", origin="lower", cmap="inferno",
                    extent=[res.z_samples[0], res.z_samples[-1], x[0], x[-1]])
    ax[0, 1].set_xlabel("z [um]")
    ax[0, 1].set_ylabel("x (lateral) [um]")
    ax[0, 1].set_title("Top view  integral|E|^2 dy")

    # (c) side view |E|^2 integrated over x  (z vs y)
    sv = res.sideview.T
    ax[1, 0].imshow(sv, aspect="auto", origin="lower", cmap="inferno",
                    extent=[res.z_samples[0], res.z_samples[-1], y[0], y[-1]])
    ax[1, 0].set_xlabel("z [um]")
    ax[1, 0].set_ylabel("y (vertical) [um]")
    ax[1, 0].set_title("Side view  integral|E|^2 dx")

    # (d) output intensity at the bus facet
    I = np.abs(res.psi_out) ** 2
    ax[1, 1].imshow(I.T, aspect="equal", origin="lower", cmap="viridis",
                    extent=[x[0], x[-1], y[0], y[-1]])
    ax[1, 1].set_xlabel("x [um]")
    ax[1, 1].set_ylabel("y [um]")
    ax[1, 1].set_title(f"Bus-facet |E|^2   (loss={ev['coupling_loss_dB']:.3f} dB, "
                       f"eff={ev['efficiency']*100:.1f}%)")

    fig.suptitle("SSC inverse taper @ 1301 nm, 2% silica, 4.5 um core", fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


# ===========================================================================
# 8. CLI demo
# ===========================================================================

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="SSC inverse-taper GA optimizer (1301nm, 2% silica)")
    ap.add_argument("--pop", type=int, default=16)
    ap.add_argument("--gen", type=int, default=12)
    ap.add_argument("--nseg", type=int, default=3)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--gds", type=str, default="ssc_best.gds")
    ap.add_argument("--fast", action="store_true", help="coarser grid, quick test")
    args = ap.parse_args()

    plat = Platform()
    print(f"Platform @ {plat.lam_um*1000:.0f} nm:  n_clad={plat.n_clad:.5f}  "
          f"n_core={plat.n_core:.5f}  t_core={plat.core_thick_um} um  "
          f"SMF w0={plat.smf_w0_um:.2f} um")

    grid = GridConfig()
    if args.fast:
        grid = GridConfig(dx_um=0.25, dy_um=0.25, dz_um=0.5,
                          x_half_um=12.0, y_half_um=10.0, save_every=12)

    bounds = GABounds(n_seg=args.nseg)
    ga = SSCGeneticOptimizer(plat, grid, bounds,
                             pop_size=args.pop, n_gen=args.gen, seed=args.seed)
    result = ga.run()

    g = result["best_gene"]
    print("\n=== BEST SSC ===")
    print(f"coupling loss @1301nm : {result['best_loss_dB']:.4f} dB")
    print(f"efficiency            : {result['best_eval']['efficiency']*100:.3f} %")
    print(f"tip width             : {g.tip_w:.3f} um")
    print(f"segment lengths [um]  : {[round(v,1) for v in g.seg_L]}")
    print(f"node widths     [um]  : {[round(v,3) for v in g.node_widths().tolist()]}")
    print(f"total length          : {g.total_length:.1f} um")

    path = gene_to_gds(g, args.gds)
    print(f"GDS written           : {path}")

    # re-evaluate best on the run grid to get full field data for plotting
    ev = evaluate_gene(g, plat, grid)
    png = plot_best(g, ev, plat, args.gds.replace(".gds", "_summary.png"))
    print(f"Plot written          : {png}")
