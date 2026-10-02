"""
ssc_ga.py -- Genetic-algorithm search of the 2.5% 3x3 um segmented SSC @ 1550 nm
==============================================================================
Fitness = insertion loss chip-mode -> SMF-28 computed by a full 3D scalar
ADI-CN BPM through the duty-averaged (effective-medium) taper, i.e. it includes
BOTH the taper radiation / non-adiabatic loss AND the facet <-> SMF mismatch.

Genes (bounds)                   DRC (min linewidth / gap = 0.5 um)
  L_um      200 .. 1500          first gap  (1-d_start)*P = 0.5  -> d_start = 1-0.5/P
  pitch_um  2.6 .. 5.0           facet tooth d_facet*P >= 0.5
  d_facet   0.12 .. 0.40         widths >= 0.5
  w_facet   3.0 .. 12.0
  p_d, p_w  0.3 .. 3.0  (power-law ramp exponents)

The effective-medium BPM is only weakly pitch-aware (through d_start and the
facet DRC clip); the pitch is therefore re-checked afterwards with the true
tooth/gap 3D BPM in verify_segmented.py.
"""
import json, math, time, sys, os
import numpy as np
from ssc_core import Platform, Grid, BPM, SSCDesign, mode_solve, db, segmented_slices

PLAT = Platform()
GRID = Grid(20.0, 18.0, 0.35)
BPMS = BPM(GRID, PLAT, pml_um=3.0)
SMF = GRID.gaussian(PLAT.smf_mfd_um / 2)
_, CHIP = mode_solve(GRID, PLAT, GRID.n2_rect(PLAT, PLAT.core_w_um))
DZ = 2.0

BOUNDS = np.array([[200, 1500], [2.6, 5.0], [0.12, 0.40], [3.0, 12.0], [0.3, 3.0], [0.3, 3.0]])
NAMES = ["L_um", "pitch_um", "d_facet", "w_facet_um", "p_d", "p_w"]


def design_from(v):
    return SSCDesign(L_um=v[0], pitch_um=v[1], d_facet=v[2], w_facet_um=v[3],
                     p_d=v[4], p_w=v[5], lead_um=20.0, W0_um=PLAT.core_w_um,
                     drc_um=PLAT.drc_min_um)


# ---- stage 2: true tooth/gap 3D BPM (pitch-aware) --------------------------
GRID2 = Grid(16.0, 14.0, float(os.environ.get("SSC_DX2", 0.15)))
BPM2 = BPM(GRID2, PLAT, pml_um=3.0)
SMF2 = GRID2.gaussian(PLAT.smf_mfd_um / 2)
_, CHIP2 = mode_solve(GRID2, PLAT, GRID2.n2_rect(PLAT, PLAT.core_w_um))
BOUNDS2 = np.array([[300, 1500], [1.15, 3.5], [0.15, 0.6], [0.8, 10.0], [0.3, 3.0], [0.3, 3.0]])
NEFF_RANGE = (PLAT.n_clad + 0.0005, 1.4623)       # local n_eff along the taper


def bragg_margin(P):
    """Distance of P / (lambda / 2 n_eff) from the nearest integer over the
    n_eff range of the taper (BPM is one-way, so contra-directional Bragg
    reflection is screened analytically instead)."""
    lo = P * 2 * NEFF_RANGE[0] / PLAT.lam_um
    hi = P * 2 * NEFF_RANGE[1] / PLAT.lam_um
    if math.floor(hi) >= math.ceil(lo):
        return 0.0
    return min(lo - math.floor(lo), math.ceil(hi) - hi)


def seg_loss(des, grid=GRID2, bpm=BPM2, chip=CHIP2, smf=SMF2, dz=0.25):
    psi, _ = bpm.run(chip.copy(), segmented_slices(des, grid, PLAT, dz))
    return db(grid.overlap_raw(psi, smf))


_cache2 = {}


def fitness2(v):
    key = tuple(np.round(v, 4))
    if key not in _cache2:
        des = design_from(v)
        pen = 0.0
        bm = bragg_margin(des.pitch_um)
        if bm < 0.08:
            pen += 2.0 * (0.08 - bm) / 0.08         # stay clear of Bragg orders
        _cache2[key] = seg_loss(des) + pen
    return _cache2[key]


def em_loss(des: SSCDesign, grid=GRID, bpm=BPMS, chip=CHIP, smf=SMF, dz=DZ):
    track = des.em_track(dz)
    slices = []
    for zl, w, d in track:
        n2 = grid.n2_rect(PLAT, w, d)
        nst = max(1, int(round(zl / dz)))
        slices.append((n2, zl / nst, nst))
    psi, _ = bpm.run(chip.copy(), slices)
    return db(grid.overlap_raw(psi, smf))


_cache = {}


def fitness(v):
    key = tuple(np.round(v, 4))
    if key not in _cache:
        _cache[key] = em_loss(design_from(v))
    return _cache[key]


def ga(pop=16, gens=14, seed=1, log=print, bounds=BOUNDS, fit=None, seeds=None):
    fit = fitness if fit is None else fit
    rng = np.random.default_rng(seed)
    lo, hi = bounds[:, 0], bounds[:, 1]
    P = lo + rng.random((pop, len(lo))) * (hi - lo)
    seeds = seeds if seeds is not None else [[700, 3.0, 0.19, 7.0, 1.0, 1.0], [1000, 3.2, 0.18, 8.0, 1.5, 1.0]]
    for i, sv in enumerate(seeds):
        P[i] = np.clip(sv, lo, hi)
    F = np.array([fit(v) for v in P])
    hist = []
    for gen in range(gens):
        order = np.argsort(F)
        P, F = P[order], F[order]
        hist.append((gen, float(F[0]), float(np.median(F))))
        log(f"gen {gen:2d}  best {F[0]:.4f} dB  median {np.median(F):.4f}  "
            + " ".join(f"{n}={x:.3g}" for n, x in zip(NAMES, P[0])))
        children = [P[0].copy(), P[1].copy()]                     # elitism
        while len(children) < pop:
            def tour():
                i, j = rng.integers(0, pop, 2)
                return P[i] if F[i] < F[j] else P[j]
            a, b = tour(), tour()
            alpha = 0.3
            cmin, cmax = np.minimum(a, b), np.maximum(a, b)
            span = cmax - cmin
            c = cmin - alpha * span + rng.random(len(lo)) * (1 + 2 * alpha) * span   # BLX-alpha
            m = rng.random(len(lo)) < 0.25
            c[m] += rng.normal(0, 0.08, m.sum()) * (hi - lo)[m]                       # mutation
            children.append(np.clip(c, lo, hi))
        P = np.array(children)
        F = np.array([fit(v) for v in P])
    order = np.argsort(F)
    return P[order], F[order], hist


if __name__ == "__main__":
    stage = sys.argv[1] if len(sys.argv) > 1 else "2"
    t0 = time.time()
    lg = lambda s: print(s, flush=True)
    if stage == "1":       # effective-medium BPM (fast, pitch-blind)
        Pp, Ff, hist = ga(log=lg)
        cache, fname = _cache, "ga_stage1_em.json"
    else:                  # true tooth/gap BPM (pitch-aware) -- the design GA
        seeds = [[800, 2.0, 0.26, 4.5, 0.7, 1.0], [1000, 2.2, 0.24, 5.0, 0.7, 1.0],
                 [1090, 2.8, 0.18, 6.5, 0.68, 1.44], [900, 1.8, 0.30, 3.5, 0.8, 1.0],
                 [1200, 2.4, 0.22, 5.5, 0.7, 1.2]]
        Pp, Ff, hist = ga(pop=16, gens=int(os.environ.get("SSC_GENS", 12)), seed=7, log=lg,
                          bounds=BOUNDS2, fit=fitness2, seeds=seeds)
        cache, fname = _cache2, "ga_stage2_seg.json"
    best = design_from(Pp[0])
    out = dict(stage=stage, best_genes=dict(zip(NAMES, map(float, Pp[0]))), best_loss_dB=float(Ff[0]),
               design=best.to_dict(), drc=best.drc(), history=hist,
               top8=[dict(zip(NAMES, map(float, v)), loss=float(f)) for v, f in zip(Pp[:8], Ff[:8])],
               evals=len(cache), runtime_s=time.time() - t0)
    json.dump(out, open(fname, "w"), indent=2)
    print(json.dumps(out["best_genes"], indent=1), out["best_loss_dB"], out["drc"])
