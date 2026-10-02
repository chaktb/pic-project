"""make_figures.py -- result figures for the 2.5% 3x3um SSC page (English labels)."""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ssc_core import Platform, Grid, mode_solve, db, d4sigma, SSCDesign

H = os.path.dirname(os.path.abspath(__file__))
V = json.load(open(os.path.join(H, "verify_results.json")))
G2 = json.load(open(os.path.join(H, "ga_stage2_seg.json")))
G3 = json.load(open(os.path.join(H, "ga_stage3_refine.json")))
G1 = json.load(open(os.path.join(H, "ga_stage1_em.json")))
F = np.load(os.path.join(H, "facet_scan.npz"))
R = np.load(os.path.join(H, "field_records.npz"))
plat = Platform()
dd = V["design"]
des = SSCDesign(L_um=dd["n_periods"] * dd["pitch_um"], pitch_um=dd["pitch_um"], d_facet=dd["d_facet"],
                w_facet_um=dd["w_facet_um"], p_d=dd["p_d"], p_w=dd["p_w"], d_start=dd["d_start"])

# ---------------- summary figure ----------------
fig = plt.figure(figsize=(13, 10))
gs = fig.add_gridspec(3, 3, hspace=0.42, wspace=0.32)

ax = fig.add_subplot(gs[0, 0])
for g, lab, c in ((G1, "Stage 1: effective-medium BPM", "#adb5bd"),
                  (G2, "Stage 2: true tooth/gap BPM", "#3b5bdb"),
                  (G3, "Stage 3: local refinement", "#e03131")):
    h = np.array(g["history"])
    off = 0 if g is not G3 else len(G2["history"])
    ax.plot(h[:, 0] + off, h[:, 1], "o-", ms=3, label=lab, color=c)
ax.set_xlabel("generation"); ax.set_ylabel("best loss (dB)"); ax.set_title("GA convergence")
ax.legend(fontsize=7); ax.grid(alpha=.3)

ax = fig.add_subplot(gs[0, 1])
Lm = np.ma.masked_invalid(F["L"])
im = ax.pcolormesh(F["D"], F["W"], Lm, shading="nearest", cmap="viridis_r", vmin=0.15, vmax=1.5)
ax.set_yscale("log"); ax.set_xlabel("facet duty"); ax.set_ylabel("facet width (um)")
ax.plot([dd["d_facet"]], [dd["w_facet_um"]], "r*", ms=12, label="GA design")
ax.axvline(0.5 / dd["pitch_um"], color="w", ls="--", lw=1)
ax.text(0.5 / dd["pitch_um"] + 0.01, 0.55, "tooth = 0.5 um\n@ P=%.2f" % dd["pitch_um"], color="w", fontsize=7)
ax.set_title("Facet-mode / SMF loss (dB)"); fig.colorbar(im, ax=ax); ax.legend(fontsize=7, loc="upper right")

ax = fig.add_subplot(gs[0, 2])
duty, width = des.periods()
z = des.lead_um + (np.arange(des.n_per) + 0.5) * des.pitch_um
ax.plot(z, duty, color="#3b5bdb", label="duty")
ax.set_xlabel("z (um)"); ax.set_ylabel("duty", color="#3b5bdb"); ax.set_ylim(0, 1.05)
ax2 = ax.twinx(); ax2.plot(z, width, color="#f08c00", label="width"); ax2.set_ylabel("width (um)", color="#f08c00")
ax2.set_ylim(0, 4)
ax.set_title("Taper law: pitch %.3f um, %d periods" % (dd["pitch_um"], dd["n_periods"]))
ax.grid(alpha=.3)

for k, (rec, coord, lab) in enumerate(((R["rx"], R["x"], "x (um), y=0"), (R["ry"], R["y"], "y (um), x=0"))):
    ax = fig.add_subplot(gs[1, k])
    A = rec / rec.max(axis=1, keepdims=True)
    ax.pcolormesh(R["z"], coord, A.T, shading="auto", cmap="magma")
    ax.set_xlabel("z (um)"); ax.set_ylabel(lab)
    ax.set_title("3D BPM |E|^2 (per-z normalized), %s cut" % ("x-z" if k == 0 else "y-z"))
    ax.set_ylim(-12, 12)

ax = fig.add_subplot(gs[1, 2])
g = Grid(16, 14, 0.10)
_, chip = mode_solve(g, plat, g.n2_rect(plat, 3.0))
smf = g.gaussian(plat.smf_mfd_um / 2)
fac = R["facet"]
ix, iy = g.nx // 2, g.ny // 2
for arr, lab, st in ((np.abs(smf), "SMF-28 (MFD 10.4 um)", "-"), (np.abs(chip), "bare 3x3 um mode", "--"),
                     (fac, "SSC output (x cut)", "-"), (fac, "SSC output (y cut)", ":")):
    cut = arr[:, iy] if "y cut" not in lab else arr[ix, :]
    co = g.x if "y cut" not in lab else g.y
    ax.plot(co, cut / cut.max(), st, label=lab)
ax.set_xlim(-12, 12); ax.set_xlabel("position (um)"); ax.set_ylabel("|E| (norm.)")
ax.set_title("Mode profiles at the facet"); ax.legend(fontsize=7); ax.grid(alpha=.3)

ax = fig.add_subplot(gs[2, 0])
wl = V["wavelength_dB"]
ax.plot([float(k) for k in wl], list(wl.values()), "o-", color="#3b5bdb")
ax.set_xlabel("wavelength (nm)"); ax.set_ylabel("coupling loss (dB)"); ax.set_title("Wavelength dependence")
ax.grid(alpha=.3)

ax = fig.add_subplot(gs[2, 1])
mis = V["misalign_dB_xy"]
dxs = [0] + [float(k) for k in mis]
ax.plot(dxs, [V["seg_bpm_loss_dB"]] + [v[0] for v in mis.values()], "o-", label="lateral (x)")
ax.plot(dxs, [V["seg_bpm_loss_dB"]] + [v[1] for v in mis.values()], "s--", label="vertical (y)")
ax.set_xlabel("SMF offset (um)"); ax.set_ylabel("coupling loss (dB)"); ax.set_title("Fiber misalignment")
ax.legend(fontsize=8); ax.grid(alpha=.3)

ax = fig.add_subplot(gs[2, 2])
cd = dict(V["cd_bias_dB"]); cd["+0.00"] = V["seg_bpm_loss_dB"]
ks = sorted(cd, key=float)
ax.plot([float(k) for k in ks], [cd[k] for k in ks], "o-", color="#e03131", label="CD bias")
cuts = V["facet_cut_dB"]
ax.plot([float(k) for k in cuts], list(cuts.values()), "s-", color="#2f9e44", label="facet cut offset")
ax.set_xlabel("CD bias / cut offset (um)"); ax.set_ylabel("coupling loss (dB)")
ax.set_title("Fabrication tolerance"); ax.legend(fontsize=8); ax.grid(alpha=.3)

fig.suptitle("2.5%% Delta 3x3 um silica -> SMF-28 @ 1550 nm, segmented SSC (DRC >= 0.5 um): "
             "%.3f dB (bare %.2f dB)" % (V["seg_bpm_loss_dB"], V["bare_loss_dB"]), fontsize=12)
fig.savefig(os.path.join(H, "ssc_result_summary.png"), dpi=120, bbox_inches="tight")
plt.close(fig)

# ---------------- method comparison bar ----------------
cal = json.load(open(os.path.join(H, "fdtd2d_calibration.json")))
fig, axs = plt.subplots(1, 2, figsize=(12, 4.2))
labs = ["bare 3x3 um\n(no SSC)", "solid inv. taper\n0.55 um tip, 2 mm", "GA seg. SSC\neff.-medium\nBPM",
        "GA seg. SSC\nlocal-mode\nEME", "GA seg. SSC\ntrue-segment\n3D BPM"]
vals = [V["bare_loss_dB"], V.get("solid_taper_dB", np.nan), V["em_bpm_loss_dB"], V["eme"]["total_dB"],
        V["seg_bpm_loss_dB"]]
cols = ["#adb5bd", "#f08c00", "#91a7ff", "#748ffc", "#3b5bdb"]
axs[0].bar(range(5), vals, color=cols)
for i, v in enumerate(vals):
    if np.isfinite(v):
        axs[0].text(i, v + 0.05, "%.3f" % v, ha="center", fontsize=9)
axs[0].set_xticks(range(5)); axs[0].set_xticklabels(labs, fontsize=8); axs[0].set_ylabel("loss (dB)")
axs[0].set_title("Coupling loss to SMF-28 @ 1550 nm")
x = np.arange(len(cal))
axs[1].bar(x - 0.18, [c["fdtd_dB"] for c in cal], 0.36, label="2D FDTD (full-wave)", color="#e03131")
axs[1].bar(x + 0.18, [c["bpm_dB"] for c in cal], 0.36, label="paraxial BPM", color="#3b5bdb")
axs[1].set_xticks(x); axs[1].set_xticklabels(["P=%.1f\nduty %.2f" % (c["pitch"], c["duty"]) for c in cal], fontsize=8)
axs[1].set_ylabel("loss of a ~60 um segmented section (dB)")
axs[1].set_title("Segment-loss calibration: 2D FDTD vs paraxial BPM"); axs[1].legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(H, "ssc_compare.png"), dpi=120); plt.close(fig)
print("figures done")
