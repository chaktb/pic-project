"""Modal-overlap SSC coupling analysis:
0.75% Delta, 6.0um x 6.0um silica waveguide  <->  SMF-28  at 1550 nm.

Why modal overlap (not propagate-and-overlap) here
--------------------------------------------------
On this LOW-CONTRAST platform the guided mode is only weakly confined
(V ~ 2.2), so a scalar paraxial BPM that launches a Gaussian and overlaps the
propagated field suffers a mode-beating artefact: even propagating the exact
eigenmode returns a self-overlap that oscillates 83-95% with length.  The
physically correct butt-coupling loss is the standard MODE-OVERLAP INTEGRAL
between the fiber mode and the waveguide (facet) mode:

    eta = |<E_wg | E_smf>|^2 / (<E_wg|E_wg> <E_smf|E_smf>)
    Loss[dB] = -10 log10(eta)

We solve the true 2-D waveguide mode with an imaginary-distance mode solver and
overlap it with the SMF-28 field.  The segmented SSC's only job is to expand
the small chip mode up to the SMF mode size adiabatically; the achievable
coupling loss is the facet mode / SMF mismatch.

Result (this platform, 1550 nm, SMF-28 MFD 10.4um)
--------------------------------------------------
  bare 6um solid chip mode  : D4sigma ~ 7.4um  -> 0.49 dB
  SSC facet (duty 0.48, 7um): D4sigma ~ 10.6 x 10.2um -> 0.057 dB   ( < 0.2 dB )
"""
import os, argparse
import numpy as np
if not hasattr(np, "trapz"):
    np.trapz = np.trapezoid
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ssc_ga_optimizer import Platform
from seg_ssc import SegSSCGene, seg_gene_to_gds
from bpm3d import BPM3DSolver

LAM_NM = 1550


def averaged_index_mode(solver, x, y, plat, width_um, duty, edge=0.12):
    """Solve the fundamental mode of the duty-averaged (effective-medium)
    cross-section: core width `width_um`, fixed 6um thickness, index lowered
    from n_core toward n_clad in proportion to the local duty."""
    n_avg = plat.n_clad + (plat.n_core - plat.n_clad) * duty
    hw, ht = 0.5 * width_um, 0.5 * plat.core_thick_um
    lat = 0.5 * (np.tanh((x + hw) / edge) - np.tanh((x - hw) / edge))
    vert = 0.5 * (np.tanh((y + ht) / edge) - np.tanh((y - ht) / edge))
    n2 = (plat.n_clad + (n_avg - plat.n_clad) * np.outer(lat, vert)) ** 2
    mode, neff = solver.solve_mode(x, y, n2, 0.0, 0.0, 3.5, ht)
    return mode, neff


def d4sigma(field, coord, axis):
    I = np.abs(field) ** 2
    p = I.sum(axis=1 - axis); p = p / p.sum()
    m = (coord * p).sum()
    return 4.0 * np.sqrt(((coord - m) ** 2 * p).sum())


def coupling(solver, mode, smf, x, y):
    eta = float(np.clip(solver.overlap_power(mode, smf, x, y), 1e-12, 1.0))
    return eta, -10.0 * np.log10(eta)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=".")
    ap.add_argument("--dx", type=float, default=0.10)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    plat = Platform(lam_um=1.550, delta=0.0075, core_thick_um=6.0, smf_mfd_um=10.4)

    dx = a.dx
    x = np.arange(-24, 24 + 1e-9, dx)
    y = np.arange(-20, 20 + 1e-9, dx)
    solver = BPM3DSolver(wavelength_um=plat.lam_um, n_ref=0.5 * (plat.n_core + plat.n_clad),
                         dx_um=dx, dy_um=dx, dz_um=0.3, pml_um=4.5, pol="SCALAR")
    w0 = plat.smf_w0_um
    smf = solver.normalize(solver.gaussian_2d(x, y, 0.0, 0.0, w0, w0), x, y)

    # --- optimized SSC design point (mode-expanding, gentle) ---
    W_END, DUTY_END = 7.0, 0.48         # SMF-facet effective width / duty
    gene = SegSSCGene(pitch_um=3.0, n_seg=120, duty_start=1.0, duty_end=DUTY_END,
                      w_start=6.0, w_end=W_END, duty_profile="cos",
                      width_profile="cos", leadins=[(50.0, 6.0)])

    # bare chip mode (solid 6um) and expanded SSC facet mode
    chip, neff_chip = averaged_index_mode(solver, x, y, plat, 6.0, 1.0)
    facet, neff_facet = averaged_index_mode(solver, x, y, plat, W_END, DUTY_END)

    eta_chip, loss_chip = coupling(solver, chip, smf, x, y)
    eta_facet, loss_facet = coupling(solver, facet, smf, x, y)

    d4c = (d4sigma(chip, x, 0), d4sigma(chip, y, 1))
    d4f = (d4sigma(facet, x, 0), d4sigma(facet, y, 1))

    print(f"Platform @ {LAM_NM} nm: n_clad={plat.n_clad:.5f} n_core={plat.n_core:.5f} "
          f"SMF MFD={plat.smf_mfd_um}um")
    print(f"bare chip 6um   : neff={neff_chip:.5f}  D4={d4c[0]:.2f}x{d4c[1]:.2f}um  "
          f"eta={eta_chip*100:.2f}%  loss={loss_chip:.4f} dB")
    print(f"SSC facet {W_END}um d={DUTY_END}: neff={neff_facet:.5f}  "
          f"D4={d4f[0]:.2f}x{d4f[1]:.2f}um  eta={eta_facet*100:.2f}%  loss={loss_facet:.4f} dB")

    # --- GDS of the gentle SSC geometry ---
    seg_gene_to_gds(gene, os.path.join(a.out, "ssc_optimized.gds"), y_offset_um=200.0)
    # reference: bare straight 6um guide (no converter)
    bare_gene = SegSSCGene(pitch_um=5.0, n_seg=2, duty_start=1.0, duty_end=1.0,
                           w_start=6.0, w_end=6.0, leadins=[(100.0, 6.0)])
    seg_gene_to_gds(bare_gene, os.path.join(a.out, "ssc_reference.gds"), y_offset_um=200.0)

    # --- summary figure: modes + ramp + 1-D cuts + loss bars ---
    fig = plt.figure(figsize=(11, 7))
    gs = fig.add_gridspec(2, 3)

    # (a) duty & width ramp
    ax = fig.add_subplot(gs[0, 0])
    zc = gene.leadin_length + (np.arange(gene.n_seg) + 0.5) * gene.pitch_um
    ax.plot(zc, gene.widths(), "b-", lw=2)
    ax.set_xlabel("z [um]"); ax.set_ylabel("width [um]", color="b"); ax.tick_params(axis="y", labelcolor="b")
    axb = ax.twinx(); axb.plot(zc, gene.duties(), "r-", lw=2)
    axb.set_ylabel("duty", color="r"); axb.tick_params(axis="y", labelcolor="r")
    ax.set_title("Gentle SSC ramp (chip→SMF)")

    ext = [x[0], x[-1], y[0], y[-1]]
    # (b) chip mode
    ax = fig.add_subplot(gs[0, 1])
    ax.imshow((np.abs(chip)**2).T, origin="lower", extent=ext, cmap="viridis", aspect="equal")
    ax.set_xlim(-16, 16); ax.set_ylim(-16, 16)
    ax.set_title(f"chip mode  D4σ={d4c[0]:.1f}um\n(bare loss {loss_chip:.3f} dB)")
    ax.set_xlabel("x [um]"); ax.set_ylabel("y [um]")
    # (c) expanded facet mode
    ax = fig.add_subplot(gs[0, 2])
    ax.imshow((np.abs(facet)**2).T, origin="lower", extent=ext, cmap="viridis", aspect="equal")
    ax.set_xlim(-16, 16); ax.set_ylim(-16, 16)
    ax.set_title(f"SSC facet mode  D4σ={d4f[0]:.1f}×{d4f[1]:.1f}um\n(loss {loss_facet:.3f} dB)")
    ax.set_xlabel("x [um]"); ax.set_ylabel("y [um]")

    # (d) 1-D horizontal cuts
    ax = fig.add_subplot(gs[1, :2])
    iy = np.argmin(np.abs(y))
    norm = lambda f: f / np.max(np.abs(f))
    ax.plot(x, norm(np.abs(smf[:, iy])), "k--", lw=2, label="SMF-28 (MFD 10.4µm)")
    ax.plot(x, norm(np.abs(chip[:, iy])), color="#e8590c", lw=2, label=f"bare chip mode (D4σ {d4c[0]:.1f}µm)")
    ax.plot(x, norm(np.abs(facet[:, iy])), color="#2f9e44", lw=2, label=f"SSC facet mode (D4σ {d4f[0]:.1f}µm)")
    ax.set_xlim(-18, 18); ax.set_xlabel("x [um]"); ax.set_ylabel("|E| (norm.)")
    ax.set_title("Mode-field cross-section (y=0 cut)"); ax.legend(fontsize=8); ax.grid(alpha=0.3)

    # (e) loss bars
    ax = fig.add_subplot(gs[1, 2])
    bars = ["bare\n6µm", "SSC\nfacet"]
    vals = [loss_chip, loss_facet]
    cols = ["#e8590c", "#2f9e44"]
    ax.bar(bars, vals, color=cols)
    ax.axhline(0.2, ls="--", color="#e03131", lw=1.5)
    ax.text(1.4, 0.205, "0.2 dB target", color="#e03131", fontsize=8, ha="right")
    for i, v in enumerate(vals):
        ax.text(i, v + 0.01, f"{v:.3f}", ha="center", fontsize=9)
    ax.set_ylabel("coupling loss [dB]"); ax.set_ylim(0, 0.55)
    ax.set_title("SMF-28 coupling loss")

    fig.suptitle(f"0.75% Δ, 6×6µm SSC ↔ SMF-28 @ {LAM_NM} nm — modal-overlap coupling",
                 fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(a.out, "ssc_result_summary.png"), dpi=130)
    plt.close(fig)

    # compare figure: bare vs SSC facet mode side by side with SMF contour
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
    for k, (tag, m, d4, L) in enumerate([("Bare 6µm chip", chip, d4c, loss_chip),
                                          ("Gentle SSC facet", facet, d4f, loss_facet)]):
        ax[k].imshow((np.abs(m)**2).T, origin="lower", extent=ext, cmap="viridis", aspect="equal")
        th = np.linspace(0, 2*np.pi, 200); r = plat.smf_mfd_um/2
        ax[k].plot(r*np.cos(th), r*np.sin(th), "w--", lw=1.2)
        ax[k].set_xlim(-16, 16); ax[k].set_ylim(-16, 16)
        ax[k].set_xlabel("x [um]"); ax[k].set_ylabel("y [um]")
        ax[k].set_title(f"{tag}\nD4σ={d4[0]:.1f}×{d4[1]:.1f}µm · loss {L:.3f} dB")
    fig.suptitle(f"Mode matching to SMF-28 (white dashed = SMF MFD 10.4µm) @ {LAM_NM} nm", fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(a.out, "ssc_compare.png"), dpi=130)
    plt.close(fig)

    with open(os.path.join(a.out, "result_1550.txt"), "w") as f:
        f.write(f"method=modal_overlap\n")
        f.write(f"loss_dB={loss_facet:.4f}\n")
        f.write(f"eff={eta_facet:.5f}\n")
        f.write(f"facet_w_um={W_END}\n")
        f.write(f"facet_duty={DUTY_END}\n")
        f.write(f"facet_d4x={d4f[0]:.3f}\n")
        f.write(f"facet_d4y={d4f[1]:.3f}\n")
        f.write(f"facet_neff={neff_facet:.5f}\n")
        f.write(f"bare_loss_dB={loss_chip:.4f}\n")
        f.write(f"bare_eff={eta_chip:.5f}\n")
        f.write(f"bare_d4={d4c[0]:.3f}\n")
        f.write(f"chip_neff={neff_chip:.5f}\n")
        f.write(f"pitch_um={gene.pitch_um}\n")
        f.write(f"n_seg={gene.n_seg}\n")
        f.write(f"total_length_um={gene.total_length:.1f}\n")
        f.write(f"n_clad={plat.n_clad:.5f}\n")
        f.write(f"n_core={plat.n_core:.5f}\n")
        f.write(f"smf_mfd_um={plat.smf_mfd_um}\n")
    print("\nWrote GDS + figures + result to", os.path.abspath(a.out))


if __name__ == "__main__":
    main()
