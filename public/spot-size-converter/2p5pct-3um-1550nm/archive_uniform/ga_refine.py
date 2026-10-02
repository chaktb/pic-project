"""Stage-3 local GA refinement around the stage-2 optimum (true segmented 3D BPM)."""
import json, time, numpy as np
import ssc_ga as G

prev = json.load(open("ga_stage2_seg.json"))
seeds = [[t[k] for k in G.NAMES] for t in prev["top8"]]
B3 = np.array([[600, 1200], [1.85, 2.35], [0.20, 0.32], [2.4, 4.5], [0.5, 2.0], [0.5, 2.0]])
t0 = time.time()
Pp, Ff, hist = G.ga(pop=12, gens=9, seed=11, log=lambda s: print(s, flush=True),
                    bounds=B3, fit=G.fitness2, seeds=seeds[:6])
best = G.design_from(Pp[0])
out = dict(stage="3-refine", best_genes=dict(zip(G.NAMES, map(float, Pp[0]))), best_loss_dB=float(Ff[0]),
           design=best.to_dict(), drc=best.drc(), history=hist,
           top8=[dict(zip(G.NAMES, map(float, v)), loss=float(f)) for v, f in zip(Pp[:8], Ff[:8])],
           evals=len(G._cache2), runtime_s=time.time() - t0)
json.dump(out, open("ga_stage3_refine.json", "w"), indent=2)
print(out["best_genes"], out["best_loss_dB"], out["drc"])
