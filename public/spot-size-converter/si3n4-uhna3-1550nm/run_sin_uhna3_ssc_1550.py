"""Si3N4 inverse-taper spot-size converter to UHNA3 fiber at 1550 nm.

Platform
--------
  3 um thermal SiO2 (BOX)  |  Si3N4 core 1.5 um (W) x 0.3 um (H)  |  5 um SiO2 top
  n_Si3N4 = 1.99628 (Luke 2015 Sellmeier),  n_SiO2 = 1.44402 (Malitson) @ 1550 nm
  Fiber: Nufern UHNA3, NA 0.35, MFD 4.1 um @ 1550 nm.

Mechanism
---------
The bright, tightly-confined 1.5 um Si3N4 mode is a poor match to the 4.1 um
fiber (~5 dB butt loss).  An INVERSE taper narrows the Si3N4 width toward a
sub-wavelength tip; as the width shrinks the effective index drops toward the
SiO2 cladding and the mode delocalises, expanding to ~fiber size.  The coupling
loss is the mode-overlap integral between the tip mode and the UHNA3 field.

Why a finite-difference eigensolver
-----------------------------------
For a high-contrast Si3N4 wire near cutoff, an imaginary-distance paraxial BPM
mode solve is unstable (it collapses onto a spurious tightly-bound mode).  We
solve the true 2-D scalar Helmholtz eigenproblem with a sparse shift-invert
eigensolve, which is robust and grid-convergent.

Result (fine grid): best tip w = 0.23 um -> mode D4sigma 4.6 x 4.5 um,
overlap 90.8% with UHNA3  ->  0.42 dB   (vs 5.0 dB for the bare 1.5 um guide).
"""
import os, math, argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from fdmode import (build_index, solve_fd, d4sigma, gauss2d, overlap,
                    N_SIN, N_OX, LAM)

W_CHIP = 1.50     # um, Si3N4 chip width
H_CORE = 0.30     # um, Si3N4 thickness
MFD = 4.1         # um, UHNA3 mode-field diameter @1550nm


def loss_of(mode, uh, x, y):
    eta = overlap(mode, uh, x, y)
    return eta, -10.0 * math.log10(max(eta, 1e-12))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=".")
    ap.add_argument("--dx", type=float, default=0.022)
    ap.add_argument("--dy", type=float, default=0.016)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    k0 = 2 * math.pi / LAM
    w0 = MFD / 2

    x = np.arange(-9, 9 + 1e-9, a.dx)
    y = np.arange(-7, 7 + 1e-9, a.dy)
    uh = gauss2d(x, y, w0)

    print(f"n_Si3N4={N_SIN:.5f}  n_SiO2={N_OX:.5f}  grid {len(x)}x{len(y)}")

    # tip-width sweep -> coupling curve
    tips = [0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.28, 0.30]
    curve = []
    best = None
    for w in tips:
        n = build_index(x, y, w)
        ne, md = solve_fd(x, y, n, k0, 1)
        m = md[0]
        eta, L = loss_of(m, uh, x, y)
        d4 = (d4sigma(m, x, 0), d4sigma(m, y, 1))
        curve.append((w, ne[0], d4[0], d4[1], eta, L))
        if best is None or L < best[-1]:
            best = (w, ne[0], d4[0], d4[1], eta, L)
        print(f"  w={w:.2f} neff={ne[0]:.5f} D4={d4[0]:.2f}x{d4[1]:.2f} "
              f"eta={eta*100:.2f}% loss={L:.4f}")
    w_opt = best[0]
    print(f"\nBEST tip w={w_opt:.2f}um  loss={best[-1]:.4f} dB (eta={best[4]*100:.2f}%)")

    # solve the optimal tip mode and the bare chip mode (high res)
    n_tip = build_index(x, y, w_opt); _, md = solve_fd(x, y, n_tip, k0, 1); tipmode = md[0]
    n_chip = build_index(x, y, W_CHIP); _, mc = solve_fd(x, y, n_chip, k0, 1); chipmode = mc[0]
    eta_tip, loss_tip = loss_of(tipmode, uh, x, y)
    eta_chip, loss_chip = loss_of(chipmode, uh, x, y)
    d4_tip = (d4sigma(tipmode, x, 0), d4sigma(tipmode, y, 1))
    d4_chip = (d4sigma(chipmode, x, 0), d4sigma(chipmode, y, 1))

    # ---- GDS: linear inverse taper on the Si3N4 layer ----
    import gdstk
    L_lead_chip, L_taper, L_tip = 20.0, 200.0, 10.0   # um
    lib = gdstk.Library(); cell = lib.new_cell("SIN_SSC")
    zc = [-L_lead_chip, 0.0, L_taper, L_taper + L_tip]
    hw = [W_CHIP/2, W_CHIP/2, w_opt/2, w_opt/2]
    top = list(zip(zc, hw)); bot = list(zip(zc[::-1], [-h for h in hw[::-1]]))
    cell.add(gdstk.Polygon(top + bot, layer=1, datatype=0))
    lib.write_gds(os.path.join(a.out, "ssc_optimized.gds"))
    # reference: straight 1.5um Si3N4 (no taper)
    lib2 = gdstk.Library(); c2 = lib2.new_cell("SIN_STRAIGHT")
    c2.add(gdstk.rectangle((-20, -W_CHIP/2), (120, W_CHIP/2), layer=1, datatype=0))
    lib2.write_gds(os.path.join(a.out, "ssc_reference.gds"))

    # ---- summary figure ----
    fig = plt.figure(figsize=(11, 7))
    gs = fig.add_gridspec(2, 3)
    ext = [x[0], x[-1], y[0], y[-1]]

    # (a) taper geometry
    ax = fig.add_subplot(gs[0, 0])
    ax.fill(np.array(zc + zc[::-1]),
            np.array(hw + [-h for h in hw[::-1]]), color="#3b5bdb", alpha=0.7)
    ax.set_xlabel("z [um]"); ax.set_ylabel("Si3N4 width [um]")
    ax.set_title(f"Inverse taper  {W_CHIP}→{w_opt:.2f} µm")
    ax.set_ylim(-1.0, 1.0); ax.grid(alpha=0.3)

    # (b) chip mode
    ax = fig.add_subplot(gs[0, 1])
    ax.imshow((np.abs(chipmode)**2).T, origin="lower", extent=ext, cmap="magma", aspect="equal")
    ax.set_xlim(-4, 4); ax.set_ylim(-4, 4)
    ax.set_title(f"chip mode {W_CHIP}µm\nD4σ={d4_chip[0]:.2f}×{d4_chip[1]:.2f}µm ({loss_chip:.2f} dB)")
    ax.set_xlabel("x [um]"); ax.set_ylabel("y [um]")

    # (c) tip mode + UHNA3 circle
    ax = fig.add_subplot(gs[0, 2])
    ax.imshow((np.abs(tipmode)**2).T, origin="lower", extent=ext, cmap="viridis", aspect="equal")
    th = np.linspace(0, 2*np.pi, 200)
    ax.plot((MFD/2)*np.cos(th), (MFD/2)*np.sin(th), "w--", lw=1.2)
    ax.set_xlim(-6, 6); ax.set_ylim(-6, 6)
    ax.set_title(f"tip mode {w_opt:.2f}µm\nD4σ={d4_tip[0]:.1f}×{d4_tip[1]:.1f}µm ({loss_tip:.3f} dB)")
    ax.set_xlabel("x [um]"); ax.set_ylabel("y [um]")

    # (d) coupling loss vs tip width
    ax = fig.add_subplot(gs[1, :2])
    ws = [c[0] for c in curve]; Ls = [c[5] for c in curve]
    ax.plot(ws, Ls, "o-", color="#3b5bdb")
    ax.axvline(w_opt, ls=":", color="#2f9e44")
    ax.annotate(f"{best[-1]:.3f} dB @ {w_opt:.2f}µm", (w_opt, best[-1]),
                textcoords="offset points", xytext=(8, 8), color="#2f9e44")
    ax.set_xlabel("Si3N4 tip width [um]"); ax.set_ylabel("UHNA3 coupling loss [dB]")
    ax.set_title("Coupling loss vs inverse-taper tip width"); ax.grid(alpha=0.3)

    # (e) 1-D mode cut
    ax = fig.add_subplot(gs[1, 2])
    iy = np.argmin(np.abs(y)); nrm = lambda f: f/np.max(np.abs(f))
    ax.plot(x, nrm(np.abs(uh[:, iy])), "k--", lw=2, label="UHNA3 (4.1µm)")
    ax.plot(x, nrm(np.abs(tipmode[:, iy])), color="#2f9e44", lw=2, label="tip mode")
    ax.plot(x, nrm(np.abs(chipmode[:, iy])), color="#e8590c", lw=1.5, label="chip mode")
    ax.set_xlim(-8, 8); ax.set_xlabel("x [um]"); ax.set_ylabel("|E| (norm.)")
    ax.legend(fontsize=8); ax.set_title("Mode-field cut (y=0)"); ax.grid(alpha=0.3)

    fig.suptitle(f"Si3N4 (1.5×0.3µm) → UHNA3 inverse-taper SSC @ {int(LAM*1000)} nm — "
                 f"modal-overlap coupling", fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(a.out, "ssc_result_summary.png"), dpi=130)
    plt.close(fig)

    # ---- compare figure ----
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
    for k, (tag, m, d4, L, cmap) in enumerate([
            (f"Bare {W_CHIP}µm chip", chipmode, d4_chip, loss_chip, "magma"),
            (f"Inverse-taper tip {w_opt:.2f}µm", tipmode, d4_tip, loss_tip, "viridis")]):
        ax[k].imshow((np.abs(m)**2).T, origin="lower", extent=ext, cmap=cmap, aspect="equal")
        ax[k].plot((MFD/2)*np.cos(th), (MFD/2)*np.sin(th), "w--", lw=1.2)
        ax[k].set_xlim(-6, 6); ax[k].set_ylim(-6, 6)
        ax[k].set_xlabel("x [um]"); ax[k].set_ylabel("y [um]")
        ax[k].set_title(f"{tag}\nD4σ={d4[0]:.2f}×{d4[1]:.2f}µm · {L:.3f} dB")
    fig.suptitle("Mode matching to UHNA3 (white dashed = 4.1µm MFD) @ 1550 nm", fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(a.out, "ssc_compare.png"), dpi=130)
    plt.close(fig)

    with open(os.path.join(a.out, "result_1550.txt"), "w") as f:
        f.write("method=fd_mode_overlap\n")
        f.write(f"loss_dB={loss_tip:.4f}\n")
        f.write(f"eff={eta_tip:.5f}\n")
        f.write(f"tip_w_um={w_opt:.3f}\n")
        f.write(f"tip_neff={best[1]:.5f}\n")
        f.write(f"tip_d4x={d4_tip[0]:.3f}\n")
        f.write(f"tip_d4y={d4_tip[1]:.3f}\n")
        f.write(f"bare_loss_dB={loss_chip:.4f}\n")
        f.write(f"bare_eff={eta_chip:.5f}\n")
        f.write(f"bare_d4x={d4_chip[0]:.3f}\n")
        f.write(f"bare_d4y={d4_chip[1]:.3f}\n")
        f.write(f"n_sin={N_SIN:.5f}\n")
        f.write(f"n_sio2={N_OX:.5f}\n")
        f.write(f"mfd_um={MFD}\n")
        f.write(f"taper_length_um={L_taper}\n")
    print("\nWrote GDS + figures + result to", os.path.abspath(a.out))


if __name__ == "__main__":
    main()
