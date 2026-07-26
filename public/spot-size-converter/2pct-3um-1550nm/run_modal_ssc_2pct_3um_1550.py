"""Modal-overlap SSC coupling analysis:
2% Delta (GeO2-doped silica), 3um x 3um core  <->  SMF-28  at 1550 nm.

Why modal overlap (not propagate-and-overlap) for the low-loss condition
------------------------------------------------------------------------
To reach <0.2 dB the chip mode must be expanded to ~9-10um (near the 9.2um SMF
mode).  On a 4.5um core that expanded facet mode is only weakly guided -- it sits
close to cutoff (n_eff only ~0.002 above the cladding).  A scalar paraxial BPM
that launches a Gaussian and overlaps the propagated field is unreliable there:
even propagating the exact expanded eigenmode returns a self-overlap of ~28%
(a mode-beating artefact), which spuriously inflates the loss to ~0.28 dB.

The physically correct butt-coupling loss is the MODE-OVERLAP INTEGRAL between
the fiber mode and the (duty-averaged effective-medium) facet mode:

    eta = |<E_wg | E_smf>|^2 / (<E_wg|E_wg> <E_smf|E_smf>)
    Loss[dB] = -10 log10(eta)

Result (this platform, 1550 nm, SMF-28 MFD 10.4um)
-------------------------------------------------
  bare 3.0um solid chip mode : D4sigma ~ 4.4um  -> ~2.6 dB
  SSC facet (duty 0.24, 7.0um): D4sigma ~11.2 x 9.8um -> ~0.15 dB  (floor; <0.1 dB
  is NOT reachable: the thin 3um core makes the expanded mode elliptical).
"""
import os, math, argparse
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
W_END, DUTY_END = 7.0, 0.24          # SMF-facet effective width / duty (low-loss point)


def averaged_index_mode(solver, x, y, plat, width_um, duty, edge=0.10):
    n_avg = plat.n_clad + (plat.n_core - plat.n_clad) * duty
    hw, ht = 0.5 * width_um, 0.5 * plat.core_thick_um
    lat = 0.5 * (np.tanh((x + hw) / edge) - np.tanh((x - hw) / edge))
    vert = 0.5 * (np.tanh((y + ht) / edge) - np.tanh((y - ht) / edge))
    n2 = (plat.n_clad + (n_avg - plat.n_clad) * np.outer(lat, vert)) ** 2
    return solver.solve_mode(x, y, n2, 0.0, 0.0, 3.5, ht)


def d4sigma(field, coord, axis):
    I = np.abs(field) ** 2
    p = I.sum(axis=1 - axis); p = p / p.sum()
    m = (coord * p).sum()
    return 4.0 * math.sqrt(((coord - m) ** 2 * p).sum())


def coupling(solver, mode, smf, x, y):
    eta = float(np.clip(solver.overlap_power(mode, smf, x, y), 1e-12, 1.0))
    return eta, -10.0 * math.log10(eta)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=".")
    ap.add_argument("--dx", type=float, default=0.09)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    plat = Platform(lam_um=1.550, delta=0.02, core_thick_um=3.0, smf_mfd_um=10.4)
    dx = a.dx
    x = np.arange(-20, 20 + 1e-9, dx)
    y = np.arange(-18, 18 + 1e-9, dx)
    solver = BPM3DSolver(wavelength_um=plat.lam_um, n_ref=0.5 * (plat.n_core + plat.n_clad),
                         dx_um=dx, dy_um=dx, dz_um=0.3, pml_um=4.5, pol="SCALAR")
    w0 = plat.smf_w0_um
    smf = solver.normalize(solver.gaussian_2d(x, y, 0.0, 0.0, w0, w0), x, y)

    # gentle mode-expanding segmented taper for the GDS
    gene = SegSSCGene(pitch_um=4.0, n_seg=120, duty_start=0.95, duty_end=DUTY_END,
                      w_start=3.0, w_end=W_END, duty_profile="cos",
                      width_profile="cos", leadins=[(50.0, 3.0)])

    chip, neff_chip = averaged_index_mode(solver, x, y, plat, 3.0, 1.0)
    facet, neff_facet = averaged_index_mode(solver, x, y, plat, W_END, DUTY_END)
    eta_chip, loss_chip = coupling(solver, chip, smf, x, y)
    eta_facet, loss_facet = coupling(solver, facet, smf, x, y)
    d4c = (d4sigma(chip, x, 0), d4sigma(chip, y, 1))
    d4f = (d4sigma(facet, x, 0), d4sigma(facet, y, 1))

    print(f"Platform @ {LAM_NM} nm: n_clad={plat.n_clad:.5f} n_core={plat.n_core:.5f} SMF MFD={plat.smf_mfd_um}um")
    print(f"bare chip 3.0um   : neff={neff_chip:.5f} D4={d4c[0]:.2f}x{d4c[1]:.2f}um eta={eta_chip*100:.2f}% loss={loss_chip:.4f} dB")
    print(f"SSC facet {W_END}um d={DUTY_END}: neff={neff_facet:.5f} (margin {neff_facet-plat.n_clad:+.5f}) "
          f"D4={d4f[0]:.2f}x{d4f[1]:.2f}um eta={eta_facet*100:.2f}% loss={loss_facet:.4f} dB")

    # duty tolerance scan around the operating point
    tol = []
    for duty in [0.18, 0.20, 0.22, 0.24, 0.26, 0.30]:
        m, ne = averaged_index_mode(solver, x, y, plat, W_END, duty)
        e, L = coupling(solver, m, smf, x, y)
        tol.append((duty, ne - plat.n_clad, e, L))
        print(f"  duty={duty:.2f} margin={ne-plat.n_clad:+.5f} eta={e*100:.2f}% loss={L:.4f} dB")

    # GDS
    seg_gene_to_gds(gene, os.path.join(a.out, "ssc_optimized.gds"), y_offset_um=200.0)
    bare_gene = SegSSCGene(pitch_um=5.0, n_seg=2, duty_start=1.0, duty_end=1.0,
                           w_start=3.0, w_end=3.0, leadins=[(100.0, 3.0)])
    seg_gene_to_gds(bare_gene, os.path.join(a.out, "ssc_reference.gds"), y_offset_um=200.0)

    # summary figure
    fig = plt.figure(figsize=(11, 7)); gs = fig.add_gridspec(2, 3)
    zc = gene.leadin_length + (np.arange(gene.n_seg) + 0.5) * gene.pitch_um
    ax = fig.add_subplot(gs[0, 0])
    ax.plot(zc, gene.widths(), "b-", lw=2); ax.set_xlabel("z [um]")
    ax.set_ylabel("width [um]", color="b"); ax.tick_params(axis="y", labelcolor="b")
    axb = ax.twinx(); axb.plot(zc, gene.duties(), "r-", lw=2)
    axb.set_ylabel("duty", color="r"); axb.tick_params(axis="y", labelcolor="r")
    ax.set_title("Mode-expanding SSC ramp (chip→SMF)")

    ext = [x[0], x[-1], y[0], y[-1]]
    ax = fig.add_subplot(gs[0, 1])
    ax.imshow((np.abs(chip)**2).T, origin="lower", extent=ext, cmap="magma", aspect="equal")
    ax.set_xlim(-8, 8); ax.set_ylim(-8, 8)
    ax.set_title(f"chip mode 3.0µm\nD4σ={d4c[0]:.1f}×{d4c[1]:.1f}µm ({loss_chip:.2f} dB)")
    ax.set_xlabel("x [um]"); ax.set_ylabel("y [um]")

    ax = fig.add_subplot(gs[0, 2])
    ax.imshow((np.abs(facet)**2).T, origin="lower", extent=ext, cmap="viridis", aspect="equal")
    th = np.linspace(0, 2*np.pi, 200)
    ax.plot((plat.smf_mfd_um/2)*np.cos(th), (plat.smf_mfd_um/2)*np.sin(th), "w--", lw=1.2)
    ax.set_xlim(-10, 10); ax.set_ylim(-10, 10)
    ax.set_title(f"SSC facet mode\nD4σ={d4f[0]:.1f}×{d4f[1]:.1f}µm ({loss_facet:.3f} dB)")
    ax.set_xlabel("x [um]"); ax.set_ylabel("y [um]")

    ax = fig.add_subplot(gs[1, :2])
    iy = np.argmin(np.abs(y)); nrm = lambda f: f/np.max(np.abs(f))
    ax.plot(x, nrm(np.abs(smf[:, iy])), "k--", lw=2, label="SMF-28 (10.4µm)")
    ax.plot(x, nrm(np.abs(facet[:, iy])), color="#2f9e44", lw=2, label=f"SSC facet (D4σ {d4f[0]:.1f}µm)")
    ax.plot(x, nrm(np.abs(chip[:, iy])), color="#e8590c", lw=1.5, label=f"bare chip (D4σ {d4c[0]:.1f}µm)")
    ax.set_xlim(-12, 12); ax.set_xlabel("x [um]"); ax.set_ylabel("|E| (norm.)")
    ax.legend(fontsize=8); ax.set_title("Mode-field cut (y=0)"); ax.grid(alpha=0.3)

    ax = fig.add_subplot(gs[1, 2])
    duties = [t[0] for t in tol]; losses = [t[3] for t in tol]
    ax.plot(duties, losses, "o-", color="#3b5bdb")
    ax.axhline(0.2, ls="--", color="#e03131", lw=1.5)
    ax.text(0.24, 0.205, "0.2 dB", color="#e03131", fontsize=8, ha="right")
    ax.set_xlabel("SMF-facet duty"); ax.set_ylabel("coupling loss [dB]")
    ax.set_title("Loss vs facet duty (w=7µm)"); ax.grid(alpha=0.3)

    fig.suptitle(f"2% Δ, 3×3µm SSC ↔ SMF-28 @ {LAM_NM} nm — modal-overlap coupling", fontsize=12)
    fig.tight_layout(); fig.savefig(os.path.join(a.out, "ssc_result_summary.png"), dpi=130); plt.close(fig)

    # compare figure
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
    for k, (tag, m, d4, L, cmap) in enumerate([
            ("Bare 3.0µm chip", chip, d4c, loss_chip, "magma"),
            ("Mode-expanding SSC facet", facet, d4f, loss_facet, "viridis")]):
        ax[k].imshow((np.abs(m)**2).T, origin="lower", extent=ext, cmap=cmap, aspect="equal")
        ax[k].plot((plat.smf_mfd_um/2)*np.cos(th), (plat.smf_mfd_um/2)*np.sin(th), "w--", lw=1.2)
        ax[k].set_xlim(-10, 10); ax[k].set_ylim(-10, 10)
        ax[k].set_xlabel("x [um]"); ax[k].set_ylabel("y [um]")
        ax[k].set_title(f"{tag}\nD4σ={d4[0]:.2f}×{d4[1]:.2f}µm · {L:.3f} dB")
    fig.suptitle("Mode matching to SMF-28 (white dashed = 10.4µm MFD) @ 1550 nm", fontsize=12)
    fig.tight_layout(); fig.savefig(os.path.join(a.out, "ssc_compare.png"), dpi=130); plt.close(fig)

    with open(os.path.join(a.out, "result_1550.txt"), "w") as f:
        f.write("method=modal_overlap\n")
        f.write(f"loss_dB={loss_facet:.4f}\n")
        f.write(f"eff={eta_facet:.5f}\n")
        f.write(f"facet_w_um={W_END}\n")
        f.write(f"facet_duty={DUTY_END}\n")
        f.write(f"facet_d4x={d4f[0]:.3f}\n")
        f.write(f"facet_d4y={d4f[1]:.3f}\n")
        f.write(f"facet_neff={neff_facet:.5f}\n")
        f.write(f"facet_margin={neff_facet-plat.n_clad:.5f}\n")
        f.write(f"bare_loss_dB={loss_chip:.4f}\n")
        f.write(f"bare_eff={eta_chip:.5f}\n")
        f.write(f"bare_d4x={d4c[0]:.3f}\n")
        f.write(f"bare_d4y={d4c[1]:.3f}\n")
        f.write(f"n_clad={plat.n_clad:.5f}\n")
        f.write(f"n_core={plat.n_core:.5f}\n")
        f.write(f"smf_mfd_um={plat.smf_mfd_um}\n")
        for duty, margin, e, L in tol:
            f.write(f"tol_duty_{duty:.2f}={L:.4f}\n")
    print("\nWrote GDS + figures + result to", os.path.abspath(a.out))


if __name__ == "__main__":
    main()
