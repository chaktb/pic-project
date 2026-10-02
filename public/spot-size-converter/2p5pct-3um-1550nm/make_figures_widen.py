"""make_figures_widen.py -- figures for the widening segmented SSC page (English labels)."""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MP
from ssc_core import Platform, Grid, mode_solve, SSCDesign
from seg_vp import from_genes, NAMES

H = os.path.dirname(os.path.abspath(__file__))
V = json.load(open(os.path.join(H, "verify_results.json")))
GA = json.load(open(os.path.join(H, "ga_widen_A.json")))
R = np.load(os.path.join(H, "field_records.npz"))
des = from_genes([V["genes"][k] for k in NAMES])
plat = Platform()
per = des.periods()
z = np.array([des.lead_um + p[0] + p[1] - p[2] / 2 for p in per])     # tooth centers
P = np.array([p[1] for p in per]); T = np.array([p[2] for p in per]); Gp = P - T
W = np.array([0.5 * (p[3] + p[4]) for p in per])
segs = des.teeth(); Ldev = segs[-1][1]

fig = plt.figure(figsize=(13, 11))
gs = fig.add_gridspec(4, 3, hspace=0.55, wspace=0.32, height_ratios=[1.0, 1, 1, 1])

ax = fig.add_subplot(gs[0, :])
for (z0, z1, w0, w1) in segs:
    ax.add_patch(MP([(z0, -w0 / 2), (z1, -w1 / 2), (z1, w1 / 2), (z0, w0 / 2)], closed=True,
                    fc="#3b5bdb", ec="none"))
ax.axvline(Ldev, color="#e03131", ls="--", lw=1)
ax.set_xlim(-5, Ldev + 10); ax.set_ylim(-4.5, 4.5); ax.set_xlabel("z (um)"); ax.set_ylabel("x (um)")
ax.set_title("GA-optimized widening segmented SSC: 3 um chip WG (left) -> %.1f um low-duty SMF facet (right), "
             "%d segments, L = %.0f um" % (des.w_facet_um, len(per), Ldev), fontsize=10)

ax = fig.add_subplot(gs[1, 0])
ax.plot(z, T / P, color="#3b5bdb", label="duty")
ax.set_ylabel("duty"); ax.set_ylim(0, 1); ax.set_xlabel("z (um)")
a2 = ax.twinx(); a2.plot(z, W, color="#f08c00", label="width"); a2.set_ylabel("segment width (um)", color="#f08c00")
a2.set_ylim(0, 8); ax.set_title("Duty and width ramps"); ax.grid(alpha=.3)

ax = fig.add_subplot(gs[1, 1])
ax.plot(z, P, color="#2f9e44", label="pitch")
ax.plot(z, T, color="#3b5bdb", lw=1, label="tooth")
ax.plot(z, Gp, color="#e03131", lw=1, label="gap")
ax.axhline(0.5, color="k", ls=":", lw=1); ax.text(z[0], 0.55, "0.5 um rule", fontsize=7)
ax.set_xlabel("z (um)"); ax.set_ylabel("length (um)"); ax.set_ylim(0, 3.5)
ax.set_title("Minimum-pitch rule: tooth & gap >= 0.5 um"); ax.legend(fontsize=7); ax.grid(alpha=.3)

ax = fig.add_subplot(gs[1, 2])
h = np.array(GA["history"])
ax.plot(h[:, 0], h[:, 1], "o-", ms=3, color="#e03131", label="best")
ax.plot(h[:, 0], h[:, 2], "s--", ms=3, color="#adb5bd", label="median")
ax.set_xlabel("generation"); ax.set_ylabel("loss (dB, dx 0.12 um)"); ax.set_ylim(0.2, 0.6)
ax.set_title("GA convergence (true-segment 3D BPM)"); ax.legend(fontsize=8); ax.grid(alpha=.3)

for k, (rec, coord, lab, ttl) in enumerate(((R["rx"], R["x"], "x (um)", "Top view |E|^2 (y = 0)"),
                                            (R["ry"], R["y"], "y (um)", "Side view |E|^2 (x = 0)"))):
    ax = fig.add_subplot(gs[2, k])
    ax.pcolormesh(R["z"], coord, rec.T, shading="auto", cmap="magma")
    ax.set_ylim(-12, 12); ax.set_xlabel("z (um)"); ax.set_ylabel(lab); ax.set_title(ttl)

ax = fig.add_subplot(gs[2, 2])
g = Grid(16, 14, 0.10)
_, chip = mode_solve(g, plat, g.n2_rect(plat, 3.0))
smf = np.exp(-(g.X ** 2 + g.Y ** 2) / 5.2 ** 2)
fac = R["facet"]; ix, iy = g.nx // 2, g.ny // 2
ax.plot(g.x, smf[:, iy], label="SMF-28 (MFD 10.4 um)")
ax.plot(g.x, np.abs(chip[:, iy]) / np.abs(chip).max(), "--", label="bare 3x3 um mode")
ax.plot(g.x, fac[:, iy] / fac.max(), label="SSC output, x cut")
ax.plot(g.y, fac[ix, :] / fac.max(), ":", label="SSC output, y cut")
ax.set_xlim(-12, 12); ax.set_xlabel("position (um)"); ax.set_ylabel("|E| (norm.)")
ax.set_title("Facet field vs SMF"); ax.legend(fontsize=7); ax.grid(alpha=.3)

ax = fig.add_subplot(gs[3, 0])
wl = V["wavelength_dB"]
ax.plot([float(k) for k in wl], list(wl.values()), "o-", color="#3b5bdb")
ax.axhline(0.3, color="#e03131", ls=":", lw=1); ax.text(1500, 0.302, "0.3 dB target", fontsize=7, color="#e03131")
ax.set_xlabel("wavelength (nm)"); ax.set_ylabel("loss (dB)"); ax.set_title("Wavelength"); ax.grid(alpha=.3)

ax = fig.add_subplot(gs[3, 1])
mis = V["misalign_dB_xy"]; d = [0] + [float(k) for k in mis]
ax.plot(d, [V["loss_dB"]] + [v[0] for v in mis.values()], "o-", label="lateral (x)")
ax.plot(d, [V["loss_dB"]] + [v[1] for v in mis.values()], "s--", label="vertical (y)")
ax.set_xlabel("SMF offset (um)"); ax.set_ylabel("loss (dB)"); ax.set_title("Fiber misalignment")
ax.legend(fontsize=8); ax.grid(alpha=.3)

ax = fig.add_subplot(gs[3, 2])
cd = dict(V["cd_bias_dB"]); cd["+0.00"] = V["loss_dB"]; ks = sorted(cd, key=float)
ax.plot([float(k) for k in ks], [cd[k] for k in ks], "o-", color="#e03131", label="CD bias")
cu = V["facet_cut_dB"]; kc = sorted(cu, key=float)
ax.plot([float(k) for k in kc], [cu[k] for k in kc], "s-", color="#2f9e44", label="facet cut offset")
ax.axhline(0.3, color="k", ls=":", lw=1)
ax.set_xlabel("CD bias / cut offset (um)"); ax.set_ylabel("loss (dB)"); ax.set_title("Fabrication tolerance")
ax.legend(fontsize=8); ax.grid(alpha=.3)

fig.suptitle("2.5%% Delta 3x3 um silica -> SMF-28 @ 1550 nm: widening segmented SSC, tooth/gap >= 0.5 um, "
             "%.3f dB (bare %.2f dB)" % (V["loss_dB"], V["bare_loss_dB"]), fontsize=12)
fig.savefig(os.path.join(H, "ssc_result_summary.png"), dpi=115, bbox_inches="tight"); plt.close(fig)

# ---- comparison: uniform-width reference vs widening optimized
old = json.load(open(os.path.join(H, "archive_uniform", "final_design_uniform.json")))["design"]
od = SSCDesign(L_um=old["n_periods"] * old["pitch_um"], pitch_um=old["pitch_um"], d_facet=old["d_facet"],
               w_facet_um=old["w_facet_um"], p_d=old["p_d"], p_w=old["p_w"], d_start=old["d_start"])
fig, axs = plt.subplots(3, 1, figsize=(12, 7), gridspec_kw=dict(height_ratios=[1, 1, 1.2]))
for ax, dd, ttl in ((axs[0], od, "Reference: uniform 3 um width, fixed pitch %.2f um -> %.3f dB" % (old["pitch_um"], V["reference_uniform_dB"])),
                    (axs[1], des, "Optimized: widening 3 -> %.1f um, minimum pitch %.2f-%.2f um -> %.3f dB"
                     % (des.w_facet_um, P.min(), P.max(), V["loss_dB"]))):
    for (z0, z1, w0, w1) in dd.teeth():
        ax.add_patch(MP([(z0, -w0 / 2), (z1, -w1 / 2), (z1, w1 / 2), (z0, w0 / 2)], closed=True, fc="#3b5bdb", ec="none"))
    ax.set_xlim(-5, 900); ax.set_ylim(-4.5, 4.5); ax.set_ylabel("x (um)"); ax.set_title(ttl, fontsize=10)
ax = axs[2]
labs = ["bare 3x3 um\n(no SSC)", "uniform-width\nsegmented (ref.)", "widening\nsegmented (GA)", "facet-mode\nlimit"]
vals = [V["bare_loss_dB"], V["reference_uniform_dB"], V["loss_dB"], V["facet_mode"]["smf_loss_dB"]]
ax.bar(range(4), vals, color=["#adb5bd", "#91a7ff", "#3b5bdb", "#ced4da"])
for i, v in enumerate(vals):
    ax.text(i, v + 0.05, "%.3f dB" % v, ha="center", fontsize=9)
ax.axhline(0.3, color="#e03131", ls=":", lw=1)
ax.set_xticks(range(4)); ax.set_xticklabels(labs, fontsize=8); ax.set_ylabel("loss (dB)"); ax.set_ylim(0, 3.7)
ax.set_title("Coupling loss to SMF-28 @ 1550 nm (true-segment 3D BPM, dx = 0.10 um)", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(H, "ssc_compare.png"), dpi=120); plt.close(fig)
print("ok")
