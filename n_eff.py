"""Effective index of a silica strip waveguide via EIM (TE-TE).

Inputs: wavelength (um), core width W (um), core height H (um),
        index contrast delta = (n_core - n_clad) / n_core * 100 (%).

n_core is computed from the Sellmeier equation (Malitson 1965) for fused silica.
"""

import math
import sys


def n_silica(wavelength_um):
    B1, B2, B3 = 0.6961663, 0.4079426, 0.8974794
    C1, C2, C3 = 0.0684043**2, 0.1162414**2, 9.896161**2
    l2 = wavelength_um**2
    n2 = 1 + B1 * l2 / (l2 - C1) + B2 * l2 / (l2 - C2) + B3 * l2 / (l2 - C3)
    return math.sqrt(n2)


def slab_neff_te(wavelength_um, thickness_um, n_core, n_clad):
    """Fundamental TE mode effective index of a symmetric slab waveguide.

    Solves u*tan(u) = w with u^2 + w^2 = V^2, V = (k0*d/2)*sqrt(n_core^2 - n_clad^2).
    """
    if n_core <= n_clad:
        return n_clad
    k0 = 2 * math.pi / wavelength_um
    V = (k0 * thickness_um / 2) * math.sqrt(n_core**2 - n_clad**2)
    if V == 0:
        return n_clad

    upper = min(math.pi / 2 - 1e-9, V - 1e-12)
    if upper <= 0:
        return n_clad

    def f(u):
        return u * math.tan(u) - math.sqrt(max(V * V - u * u, 0.0))

    lo, hi = 1e-9, upper
    flo, fhi = f(lo), f(hi)
    if flo * fhi > 0:
        return n_clad

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
    kx = 2 * u / thickness_um
    beta2 = (n_core * k0) ** 2 - kx * kx
    return math.sqrt(beta2) / k0


def n_eff(wavelength_um, width_um, height_um, delta_percent):
    n_core = n_silica(wavelength_um)
    n_clad = n_core * (1 - delta_percent / 100.0)
    n_eff_v = slab_neff_te(wavelength_um, height_um, n_core, n_clad)
    n_eff_total = slab_neff_te(wavelength_um, width_um, n_eff_v, n_clad)
    return {
        "n_core": n_core,
        "n_clad": n_clad,
        "n_eff_vertical": n_eff_v,
        "n_eff": n_eff_total,
    }


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python n_eff.py <wavelength_um> <width_um> <height_um> <delta_percent>")
        print("Example: python n_eff.py 1.55 0.5 0.22 0.5")
        sys.exit(1)
    wl, w, h, d = (float(x) for x in sys.argv[1:5])
    r = n_eff(wl, w, h, d)
    print(f"wavelength            = {wl:.4f} um")
    print(f"core size  W x H      = {w:.4f} x {h:.4f} um")
    print(f"index contrast        = {d:.4f} %")
    print(f"n_core (Sellmeier)    = {r['n_core']:.6f}")
    print(f"n_clad                = {r['n_clad']:.6f}")
    print(f"n_eff (vertical slab) = {r['n_eff_vertical']:.6f}")
    print(f"n_eff (full EIM)      = {r['n_eff']:.6f}")
