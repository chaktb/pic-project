"""Facet-mode landscape: SMF-28 overlap loss vs facet width and duty
(effective-medium cross-section, scalar FD mode solve)."""
import numpy as np, json
from ssc_core import *
p = Platform(); g = Grid(22, 20, 0.25); smf = g.gaussian(p.smf_mfd_um / 2)
W = np.array([0.5,0.6,0.7,0.8,1.0,1.25,1.5,2,2.5,3,3.5,4,5,6,7,8,9,10,12])
D = np.array([0.12,0.15,0.18,0.2,0.22,0.25,0.28,0.32,0.36,0.4,0.45,0.5,0.6,0.7,0.8,0.9,1.0])
L = np.full((len(W), len(D)), np.nan); M = np.full_like(L, np.nan)
for i, w in enumerate(W):
    for j, d in enumerate(D):
        if w * d < 0.3 * 0.5 or (d < 1 and False): pass
        ne, E = mode_solve(g, p, g.n2_rect(p, w, d))
        M[i, j] = ne - p.n_clad
        if M[i, j] > 2e-4:
            L[i, j] = db(g.overlap(E, smf))
k = np.nanargmin(L); i, j = np.unravel_index(k, L.shape)
print("best: w=%.2f d=%.2f loss=%.4f dB margin=%.5f" % (W[i], D[j], L[i, j], M[i, j]))
for i, w in enumerate(W):
    print("%5.2f " % w + " ".join("%5.2f" % v if np.isfinite(v) else "  -- " for v in L[i]))
np.savez("facet_scan.npz", W=W, D=D, L=L, M=M)
