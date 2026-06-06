"""WDM directional coupler — spectral response over a wavelength range.

For a fixed-geometry DC (W, H, Δ, gap, L_DC), kappa(lambda) varies with
wavelength because n_core (Sellmeier), n_eff and gamma all depend on lambda.
So P_cross(lambda) = sin^2(kappa(lambda) * L_DC) oscillates with lambda,
which is the WDM filtering action.

Usage: python wdm_coupler.py <wl_start_um> <wl_end_um> <W_um> <H_um>
                              <delta_pct> <gap_um> <L_DC_um> [N_pts=200]
"""

import math
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def n_silica(wl):
    B1, B2, B3 = 0.6961663, 0.4079426, 0.8974794
    C1, C2, C3 = 0.0684043**2, 0.1162414**2, 9.896161**2
    l2 = wl * wl
    return math.sqrt(1 + B1*l2/(l2-C1) + B2*l2/(l2-C2) + B3*l2/(l2-C3))


def slab_te_fundamental(wl, d, n_core, n_clad):
    if n_core <= n_clad:
        return None
    k0 = 2 * math.pi / wl
    a = d / 2
    V = k0 * a * math.sqrt(n_core**2 - n_clad**2)
    upper = min(math.pi/2 - 1e-9, V - 1e-12)
    if upper <= 0:
        return None

    def f(u):
        w = math.sqrt(max(V*V - u*u, 0.0))
        return u * math.tan(u) - w

    lo, hi = 1e-9, upper
    flo, fhi = f(lo), f(hi)
    if flo * fhi > 0:
        return None
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        fmid = f(mid)
        if flo * fmid <= 0:
            hi, fhi = mid, fmid
        else:
            lo, flo = mid, fmid
        if hi - lo < 1e-13:
            break
    u = 0.5 * (lo + hi)
    w_v = math.sqrt(max(V*V - u*u, 0.0))
    return u/a, w_v/a, math.sqrt(n_core**2 - (u/a/k0)**2)


def kappa_at(wl, W, H, delta, gap):
    n_clad = n_silica(wl)
    n_core = n_clad / (1 - delta / 100)
    vmode = slab_te_fundamental(wl, H, n_core, n_clad)
    if vmode is None:
        return None
    _, _, n_eff_v = vmode
    hmode = slab_te_fundamental(wl, W, n_eff_v, n_clad)
    if hmode is None:
        return None
    kappa_x, gamma_x, n_eff = hmode
    k0 = 2 * math.pi / wl
    beta = n_eff * k0
    a = W / 2
    h2 = kappa_x * kappa_x
    return 2 * h2 * gamma_x * math.exp(-gamma_x * gap) / (beta * (h2 + gamma_x*gamma_x) * (1 + gamma_x * a))


def main():
    if len(sys.argv) not in (8, 9):
        print("Usage: python wdm_coupler.py <wl_start> <wl_end> <W> <H> <delta_pct> <gap> <L_DC> [N=200]")
        print("Example: python wdm_coupler.py 1.50 1.60 5.0 5.0 0.75 5.0 600")
        sys.exit(1)
    wl_a = float(sys.argv[1])
    wl_b = float(sys.argv[2])
    W    = float(sys.argv[3])
    H    = float(sys.argv[4])
    dl   = float(sys.argv[5])
    gap  = float(sys.argv[6])
    L_DC = float(sys.argv[7])
    N    = int(sys.argv[8]) if len(sys.argv) >= 9 else 200

    wls = np.linspace(wl_a, wl_b, N)
    Pc, Pt = [], []
    for wl in wls:
        k = kappa_at(float(wl), W, H, dl, gap)
        if k is None:
            Pc.append(np.nan); Pt.append(np.nan); continue
        Pc.append(math.sin(k * L_DC) ** 2)
        Pt.append(math.cos(k * L_DC) ** 2)

    print(f"# WDM DC spectrum  W={W} H={H} delta={dl}% gap={gap} L_DC={L_DC}")
    print("# wavelength_um, P_cross, P_through")
    for wl, pc, pt in zip(wls, Pc, Pt):
        print(f"{wl:.6f}, {pc:.6f}, {pt:.6f}")

    plt.figure(figsize=(8, 4.5))
    plt.plot(wls, Pt, color="tab:blue", lw=2, label="P_through")
    plt.plot(wls, Pc, color="tab:red",  lw=2, label="P_cross")
    plt.axhline(0.5, color="gray", lw=0.7, ls="--")
    plt.xlabel("Wavelength (μm)")
    plt.ylabel("Power")
    plt.title(f"WDM DC spectrum  (W={W}μm H={H}μm Δ={dl}% gap={gap}μm L={L_DC}μm)")
    plt.ylim(0, 1.05)
    plt.grid(True, ls="--", alpha=0.5)
    plt.legend(loc="upper right")
    plt.tight_layout()
    plt.savefig("wdm_spectrum.png", dpi=150)
    sys.stderr.write("Saved: wdm_spectrum.png\n")


if __name__ == "__main__":
    main()
