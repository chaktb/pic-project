"""GA + 3D-BPM search for the lowest-coupling-loss segmented SSC
for a 2% Delta, 4.5um x 4.5um silica waveguide coupling to SMF-28 at 1550 nm.

Chip facet fixed at 4.5um solid; SMF facet fixed at 8.5um low-duty.
GA freely tunes pitch, segment count, duty ramp and width ramp profile.

Outputs (into --out dir):
  ssc_optimized.gds     GA-best mask layout
  ssc_reference.gds     un-optimized reference segmented SSC
  ssc_result_summary.png  geometry / ramps / top+side BPM / chip-facet mode
  ssc_compare.png       reference vs optimized field + summary
"""
import os, sys, time, argparse
import numpy as np
if not hasattr(np, "trapz"):
    np.trapz = np.trapezoid          # numpy>=2.0 shim (bpm3d uses np.trapz)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ssc_ga_optimizer import Platform, GridConfig
from seg_ssc import (SegSSCGene, SegGABounds, SegSSCGeneticOptimizer,
                     evaluate_seg_gene, seg_gene_to_gds)

LAM_NM = 1550


def plot_summary(gene, ev, plat, path):
    res = ev["result"]
    x, y = ev["x"], ev["y"]
    fig, ax = plt.subplots(2, 2, figsize=(11, 7))

    # (a) duty + width ramp along the dashed region
    t = gene.seg_index_t()
    zc = gene.leadin_length + (np.arange(gene.n_seg) + 0.5) * gene.pitch_um
    a0 = ax[0, 0]
    a0.plot(zc, gene.widths(), "b-o", ms=3, label="segment width [um]")
    a0.set_xlabel("z (propagation) [um]"); a0.set_ylabel("width [um]", color="b")
    a0.tick_params(axis="y", labelcolor="b")
    a0b = a0.twinx()
    a0b.plot(zc, gene.duties(), "r-s", ms=3, label="duty")
    a0b.set_ylabel("duty (fill/pitch)", color="r"); a0b.tick_params(axis="y", labelcolor="r")
    a0.set_title(f"Duty & width ramp  (pitch={gene.pitch_um:.2f}um, n={gene.n_seg})")
    a0.grid(alpha=0.3)

    # (b) top view |E|^2 integrated over y  (z vs x)
    tv = res.topview.T
    ax[0, 1].imshow(tv, aspect="auto", origin="lower", cmap="inferno",
                    extent=[res.z_samples[0], res.z_samples[-1], x[0], x[-1]])
    ax[0, 1].set_xlabel("z [um]"); ax[0, 1].set_ylabel("x (lateral) [um]")
    ax[0, 1].set_title("Top view  ∫|E|² dy")

    # (c) side view |E|^2 integrated over x  (z vs y)
    sv = res.sideview.T
    ax[1, 0].imshow(sv, aspect="auto", origin="lower", cmap="inferno",
                    extent=[res.z_samples[0], res.z_samples[-1], y[0], y[-1]])
    ax[1, 0].set_xlabel("z [um]"); ax[1, 0].set_ylabel("y (vertical) [um]")
    ax[1, 0].set_title("Side view  ∫|E|² dx")

    # (d) chip-facet output intensity
    I = np.abs(res.psi_out) ** 2
    ax[1, 1].imshow(I.T, aspect="equal", origin="lower", cmap="viridis",
                    extent=[x[0], x[-1], y[0], y[-1]])
    ax[1, 1].set_xlabel("x [um]"); ax[1, 1].set_ylabel("y [um]")
    ax[1, 1].set_title(f"Chip-facet |E|²  (loss={ev['coupling_loss_dB']:.3f} dB, "
                       f"eff={ev['efficiency']*100:.1f}%)")

    fig.suptitle(f"Segmented SSC @ {LAM_NM} nm  ·  2% Δ silica  ·  4.5 × 4.5 um core",
                 fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


def plot_compare(ref_gene, ref_ev, opt_gene, opt_ev, plat, path):
    fig, ax = plt.subplots(2, 2, figsize=(11, 7))
    for row, (tag, g, ev) in enumerate(
            [("Reference", ref_gene, ref_ev), ("GA-optimized", opt_gene, opt_ev)]):
        res = ev["result"]; x, y = ev["x"], ev["y"]
        tv = res.topview.T
        ax[row, 0].imshow(tv, aspect="auto", origin="lower", cmap="inferno",
                          extent=[res.z_samples[0], res.z_samples[-1], x[0], x[-1]])
        ax[row, 0].set_ylabel("x [um]")
        ax[row, 0].set_title(f"{tag}: top view  "
                             f"(loss={ev['coupling_loss_dB']:.3f} dB, "
                             f"pitch={g.pitch_um:.2f}um, n={g.n_seg}, L={g.total_length:.0f}um)")
        I = np.abs(res.psi_out) ** 2
        ax[row, 1].imshow(I.T, aspect="equal", origin="lower", cmap="viridis",
                          extent=[x[0], x[-1], y[0], y[-1]])
        ax[row, 1].set_ylabel("y [um]")
        ax[row, 1].set_title(f"{tag}: chip-facet |E|²  (eff={ev['efficiency']*100:.1f}%)")
    ax[1, 0].set_xlabel("z [um]"); ax[1, 1].set_xlabel("x [um]")
    fig.suptitle(f"Segmented SSC @ {LAM_NM} nm: reference vs GA-optimized "
                 f"(2% Δ, 4.5um core)", fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pop", type=int, default=18)
    ap.add_argument("--gen", type=int, default=14)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default=".")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    # 2% Delta silica @1550nm; SMF-28 MFD ~10.4um at 1550nm
    plat = Platform(lam_um=1.550, delta=0.02, core_thick_um=4.5, smf_mfd_um=10.4)

    # GA search grid (fast) and validation grid (fine)
    grid = GridConfig(dx_um=0.25, dy_um=0.25, dz_um=0.6,
                      x_half_um=16, y_half_um=13, pml_um=3.0, save_every=20)
    fine = GridConfig(dx_um=0.15, dy_um=0.15, dz_um=0.35,
                      x_half_um=18, y_half_um=15, pml_um=3.5, save_every=12)

    bounds = SegGABounds(
        pitch_um=(3.0, 8.0),
        duty_start=(0.75, 0.98),
        duty_end=(0.10, 0.40),
        w_start=(4.5, 4.5),         # chip facet fixed solid 4.5um
        w_end=8.5,                  # SMF facet fixed 8.5um
        n_seg=(35, 100),
        profiles=("linear", "quad", "cos"),
    )
    print(f"Platform @ {LAM_NM} nm: n_clad={plat.n_clad:.5f} n_core={plat.n_core:.5f} "
          f"t={plat.core_thick_um}um SMF MFD={plat.smf_mfd_um}um (w0={plat.smf_w0_um:.2f}um)")
    print("Fixed: chip=4.5um solid,  SMF facet=8.5um,  thickness=4.5um\n")

    ga = SegSSCGeneticOptimizer(plat, grid, bounds, pop_size=a.pop, n_gen=a.gen,
                                seed=a.seed, leadins=[(50.0, 4.5)])
    t0 = time.time()
    resGA = ga.run()
    g = resGA["best_gene"]
    print(f"\nGA done in {time.time()-t0:.0f}s")

    # re-validate the winner on the fine grid
    print("Re-evaluating best on fine grid ...")
    opt_ev = evaluate_seg_gene(g, plat, fine)
    print("\n=== BEST SEGMENTED SSC @ 1550 nm ===")
    print(f"coupling loss (fine)  : {opt_ev['coupling_loss_dB']:.4f} dB "
          f"(eff={opt_ev['efficiency']*100:.3f} %)")
    print(f"pitch / n_seg         : {g.pitch_um:.3f} um / {g.n_seg}")
    print(f"duty  start->end      : {g.duty_start:.3f} -> {g.duty_end:.3f} ({g.duty_profile})")
    print(f"width start->end      : {g.w_start:.3f} -> {g.w_end:.3f} um ({g.width_profile})")
    print(f"total length          : {g.total_length:.1f} um")
    print(f"residual power        : {opt_ev['residual_power']:.3f}")

    # reference (un-optimized) segmented SSC for comparison
    ref = SegSSCGene(pitch_um=5.19, n_seg=60, duty_start=0.79, duty_end=0.25,
                     w_start=4.5, w_end=8.5, duty_profile="linear",
                     width_profile="linear", leadins=[(50.0, 4.5)])
    ref_ev = evaluate_seg_gene(ref, plat, fine)
    print(f"\nReference loss (fine) : {ref_ev['coupling_loss_dB']:.4f} dB "
          f"(eff={ref_ev['efficiency']*100:.2f} %)")

    # write GDS + figures
    seg_gene_to_gds(g, os.path.join(a.out, "ssc_optimized.gds"), y_offset_um=200.0)
    seg_gene_to_gds(ref, os.path.join(a.out, "ssc_reference.gds"), y_offset_um=200.0)
    plot_summary(g, opt_ev, plat, os.path.join(a.out, "ssc_result_summary.png"))
    plot_compare(ref, ref_ev, g, opt_ev, plat, os.path.join(a.out, "ssc_compare.png"))
    print("\nWrote GDS + figures to", os.path.abspath(a.out))

    # dump a small machine-readable result for the HTML
    with open(os.path.join(a.out, "result_1550.txt"), "w") as f:
        f.write(f"loss_dB={opt_ev['coupling_loss_dB']:.4f}\n")
        f.write(f"eff={opt_ev['efficiency']:.5f}\n")
        f.write(f"pitch_um={g.pitch_um:.4f}\n")
        f.write(f"n_seg={g.n_seg}\n")
        f.write(f"duty_start={g.duty_start:.4f}\n")
        f.write(f"duty_end={g.duty_end:.4f}\n")
        f.write(f"duty_profile={g.duty_profile}\n")
        f.write(f"width_profile={g.width_profile}\n")
        f.write(f"w_start={g.w_start:.4f}\n")
        f.write(f"w_end={g.w_end:.4f}\n")
        f.write(f"total_length_um={g.total_length:.2f}\n")
        f.write(f"residual={opt_ev['residual_power']:.4f}\n")
        f.write(f"ref_loss_dB={ref_ev['coupling_loss_dB']:.4f}\n")
        f.write(f"ref_eff={ref_ev['efficiency']:.5f}\n")
        f.write(f"n_clad={plat.n_clad:.5f}\n")
        f.write(f"n_core={plat.n_core:.5f}\n")
        f.write(f"smf_mfd_um={plat.smf_mfd_um}\n")


if __name__ == "__main__":
    main()
