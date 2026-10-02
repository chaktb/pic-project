"""Reference alternative: solid lateral inverse taper (no segmentation), 3 um -> narrow tip,
DRC tip >= 0.5 um, same fine-grid BPM.  Best of a small tip/length scan is stored."""
import json, numpy as np
from ssc_core import *
p = Platform(); g = Grid(16, 14, 0.10); b = BPM(g, p); smf = g.gaussian(5.2)
_, chip = mode_solve(g, p, g.n2_rect(p, 3.0))
res = []
for wt in (0.5, 0.55, 0.6, 0.65):
    for L in (1200, 2000):
        n = int(L)
        sl = [(g.n2_rect(p, 3.0), 1.0, 20)]
        for k in range(n):
            s = (k + 0.5) / n
            w = 3.0 + (wt - 3.0) * s ** 0.5      # fast start, slow near the tip
            sl.append((g.n2_rect(p, w), 1.0, 1))
        psi, _ = b.run(chip.copy(), sl)
        res.append(dict(tip_um=wt, L_um=L, loss_dB=db(g.overlap_raw(psi, smf))))
        print(res[-1], flush=True)
best = min(res, key=lambda r: r["loss_dB"])
V = json.load(open("verify_results.json")); V["solid_taper_dB"] = best["loss_dB"]; V["solid_taper_best"] = best
V["solid_taper_scan"] = res
json.dump(V, open("verify_results.json", "w"), indent=2, default=float)
print("best", best)
