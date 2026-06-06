"""V-b curve for a symmetric step-index slab waveguide.

V = (k0 * a) * sqrt(n_core^2 - n_clad^2),  a = d/2  (half-thickness)
b = (n_eff^2 - n_clad^2) / (n_core^2 - n_clad^2)

TE_m cutoff at V_c = m * pi/2.  Single-mode condition: V < pi/2 (only TE_0).

Usage:  python vb_curve.py <wavelength_um> <thickness_um> <delta_percent>
"""

import math
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def n_silica(wavelength_um):
    B1, B2, B3 = 0.6961663, 0.4079426, 0.8974794
    C1, C2, C3 = 0.0684043**2, 0.1162414**2, 9.896161**2
    l2 = wavelength_um**2
    n2 = 1 + B1 * l2 / (l2 - C1) + B2 * l2 / (l2 - C2) + B3 * l2 / (l2 - C3)
    return math.sqrt(n2)


def b_for_mode(V, m):
    """Normalized propagation constant b for TE_m at given V (None if below cutoff)."""
    Vc = m * math.pi / 2
    if V <= Vc:
        return None
    is_even = (m % 2 == 0)
    u_min = max(Vc + 1e-9, 1e-9)
    u_max = min((m + 1) * math.pi / 2 - 1e-9, V - 1e-12)
    if u_max <= u_min:
        return None

    def f(u):
        w = math.sqrt(max(V * V - u * u, 0.0))
        if is_even:
            return u * math.tan(u) - w
        else:
            return -u / math.tan(u) - w

    lo, hi = u_min, u_max
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
        if hi - lo < 1e-12:
            break
    u = 0.5 * (lo + hi)
    w = math.sqrt(max(V * V - u * u, 0.0))
    return (w / V) ** 2


def operating_V(wavelength_um, thickness_um, n_core, n_clad):
    if n_core <= n_clad:
        return 0.0
    k0 = 2 * math.pi / wavelength_um
    a = thickness_um / 2
    return k0 * a * math.sqrt(n_core**2 - n_clad**2)


def main():
    if len(sys.argv) != 4:
        print("Usage: python vb_curve.py <wavelength_um> <thickness_um> <delta_percent>")
        print("Example: python vb_curve.py 1.55 5.0 0.75")
        sys.exit(1)
    wl, d, delta = (float(x) for x in sys.argv[1:4])

    n_clad = n_silica(wl)
    n_core = n_clad / (1 - delta / 100)
    V_op = operating_V(wl, d, n_core, n_clad)
    single_mode = V_op < math.pi / 2

    print(f"wavelength      = {wl:.4f} um")
    print(f"core thickness  = {d:.4f} um")
    print(f"index contrast  = {delta:.4f} %")
    print(f"n_core (Sellmeier) = {n_core:.6f}")
    print(f"n_clad             = {n_clad:.6f}")
    print(f"V_operating        = {V_op:.4f}")
    print(f"V_cutoff (TE_1)    = {math.pi/2:.4f}")
    print(f"single mode        = {'YES' if single_mode else 'NO'}")

    # Plot
    V_max = max(10.0, V_op + 2)
    Vs = np.linspace(0.001, V_max, 600)
    fig, ax = plt.subplots(figsize=(8, 5))

    ax.axvspan(0, math.pi / 2, color="#e3f2fd", alpha=0.6,
               label="Single-mode region (V < π/2)")

    colors = ["tab:red", "tab:blue", "tab:green", "tab:orange"]
    for m in range(4):
        bs = [b_for_mode(V, m) for V in Vs]
        Vs_plot = [V for V, b in zip(Vs, bs) if b is not None]
        bs_plot = [b for b in bs if b is not None]
        if bs_plot:
            ax.plot(Vs_plot, bs_plot, color=colors[m], linewidth=2, label=f"TE_{m}")

    # Operating point
    ax.axvline(V_op, color="black", linestyle="--", linewidth=1.5,
               label=f"V_op = {V_op:.2f}")
    b0_op = b_for_mode(V_op, 0)
    if b0_op is not None:
        ax.plot([V_op], [b0_op], "ko", markersize=8)
        ax.annotate(f"  b₀ = {b0_op:.3f}", (V_op, b0_op))

    ax.set_xlim(0, V_max)
    ax.set_ylim(0, 1)
    ax.set_xlabel("Normalized frequency V")
    ax.set_ylabel("Normalized propagation constant b")
    ax.set_title(f"V-b curve  (λ={wl} μm, d={d} μm, Δ={delta}%) "
                 f"— {'single mode' if single_mode else 'multi mode'}")
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig("vb_curve.png", dpi=150)
    print("Saved: vb_curve.png")


if __name__ == "__main__":
    main()
