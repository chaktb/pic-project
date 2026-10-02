"""
ga_widen.py -- GA for the WIDENING, minimum-pitch segmented SSC (seg_vp.WideningSSC)
Fitness = chip-mode -> SMF-28 loss from the scalar ADI-CN 3D BPM through the real
tooth/gap geometry (all teeth and gaps >= 0.5 um by construction).
Genes: L, d_start, d_facet, W_facet (>= 4.5 um: forced widening), p_d, p_w, P_min
"""
import json, os, sys, time
import numpy as np
from ssc_core import Platform, Grid, BPM, mode_solve, db, segmented_slices
from seg_vp import from_genes, NAMES

PLAT = Platform()
GRID = Grid(16.0, 14.0, float(os.environ.get("SSC_DX", 0.12)))
BPMS = BPM(GRID, PLAT, pml_um=3.0)
SMF = GRID.gaussian(PLAT.smf_mfd_um / 2)
_, CHIP = mode_solve(GRID, PLAT, GRID.n2_rect(PLAT, PLAT.core_w_um))
BOUNDS = np.array([[400, 1500], [0.70, 0.90], [0.15, 0.35], [4.5, 9.0], [0.4, 2.5], [0.4, 2.5], [1.15, 1.5]])
_cache = {}


def loss(v, grid=GRID, bpm=BPMS, chip=CHIP, smf=SMF, dz=0.25):
    psi, _ = bpm.run(chip.copy(), segmented_slices(from_genes(v), grid, PLAT, dz))
    return db(grid.overlap_raw(psi, smf))


def fitness(v):
    k = tuple(np.round(v, 4))
    if k not in _cache:
        _cache[k] = loss(v)
    return _cache[k]


def ga(pop=16, gens=14, seed=3, seeds=(), bounds=BOUNDS, log=print):
    rng = np.random.default_rng(seed)
    lo, hi = bounds[:, 0], bounds[:, 1]
    P = lo + rng.random((pop, len(lo))) * (hi - lo)
    for i, s in enumerate(seeds):
        P[i] = np.clip(s, lo, hi)
    F = np.array([fitness(v) for v in P]); hist = []
    for gen in range(gens):
        o = np.argsort(F); P, F = P[o], F[o]
        hist.append((gen, float(F[0]), float(np.median(F))))
        log(f"gen {gen:2d} best {F[0]:.4f} dB median {np.median(F):.4f} " +
            " ".join(f"{n}={x:.3g}" for n, x in zip(NAMES, P[0])))
        kids = [P[0].copy(), P[1].copy()]
        while len(kids) < pop:
            def tour():
                i, j = rng.integers(0, pop, 2); return P[i] if F[i] < F[j] else P[j]
            a, b = tour(), tour()
            mn, mx = np.minimum(a, b), np.maximum(a, b); sp = mx - mn
            c = mn - 0.3 * sp + rng.random(len(lo)) * 1.6 * sp
            m = rng.random(len(lo)) < 0.25
            c[m] += rng.normal(0, 0.08, m.sum()) * (hi - lo)[m]
            kids.append(np.clip(c, lo, hi))
        P = np.array(kids); F = np.array([fitness(v) for v in P])
    o = np.argsort(F)
    return P[o], F[o], hist


if __name__ == "__main__":
    stage = sys.argv[1] if len(sys.argv) > 1 else "A"
    t0 = time.time(); lg = lambda s: print(s, flush=True)
    if stage == "A":
        seeds = [[800, 0.8, 0.22, 6, 1, 1, 1.2], [600, 0.8, 0.22, 7, 1, 1, 1.2],
                 [1000, 0.8, 0.2, 7, 0.8, 1, 1.2], [700, 0.75, 0.24, 6, 1, 0.7, 1.25]]
        Pp, Ff, hist = ga(seeds=seeds, log=lg); fn = "ga_widen_A.json"
    else:
        prev = json.load(open("ga_widen_A.json"))
        seeds = [[t[k] for k in NAMES] for t in prev["top8"]][:6]
        b = np.array(seeds); lo = b.min(0); hi = b.max(0); span = np.maximum(hi - lo, 0.05 * (BOUNDS[:, 1] - BOUNDS[:, 0]))
        B = np.stack([np.maximum(BOUNDS[:, 0], lo - 0.5 * span), np.minimum(BOUNDS[:, 1], hi + 0.5 * span)], 1)
        Pp, Ff, hist = ga(pop=12, gens=10, seed=5, seeds=seeds, bounds=B, log=lg); fn = "ga_widen_B.json"
    best = from_genes(Pp[0])
    out = dict(stage=stage, best_genes=dict(zip(NAMES, map(float, Pp[0]))), best_loss_dB=float(Ff[0]),
               design=best.to_dict(), drc=best.drc(), history=hist, grid_dx=GRID.dx,
               top8=[dict(zip(NAMES, map(float, v)), loss=float(f)) for v, f in zip(Pp[:8], Ff[:8])],
               evals=len(_cache), runtime_s=time.time() - t0)
    json.dump(out, open(fn, "w"), indent=2)
    print(out["best_genes"], out["best_loss_dB"], out["drc"])
