"""Directional coupler (two identical parallel waveguides) design.

Coupled-mode coupling coefficient (Yariv form, slab-style):

    kappa = 2 * h^2 * q * exp(-q * gap) / ( beta * (h^2 + q^2) * (1 + q * a) )

with h = kappa_x (transverse k inside core), q = gamma_x (cladding decay),
beta = n_eff * k0, a = W/2, all from the horizontal slab fundamental TE mode
(EIM: vertical first, horizontal second).

Power transfer (lossless):
    P_cross(z)   = sin^2(kappa * z)
    P_through(z) = cos^2(kappa * z)
    L_c          = pi / (2 * kappa)        (full transfer)
    L_50         = pi / (4 * kappa)        (3 dB / 50:50)

Loss(dB) = -10 log10(power fraction) at each output port.
Sellmeier (Malitson 1965) gives n_core(lambda) for fused silica.
"""

import math
import sys


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


def coupling_coefficient(wl, W, H, delta, gap):
    n_core = n_silica(wl)
    n_clad = n_core * (1 - delta / 100)
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
    q = gamma_x
    kappa = 2 * h2 * q * math.exp(-q * gap) / (beta * (h2 + q*q) * (1 + q * a))
    return {
        "n_core": n_core, "n_clad": n_clad, "n_eff": n_eff,
        "kappa_x": kappa_x, "gamma_x": gamma_x, "beta": beta,
        "kappa": kappa,
    }


def main():
    if len(sys.argv) != 7:
        print("Usage: python directional_coupler.py <wl_um> <W_um> <H_um> <delta_pct> <gap_um> <L_DC_um>")
        print("Example: python directional_coupler.py 1.55 5.0 5.0 0.75 3.0 200")
        sys.exit(1)
    wl    = float(sys.argv[1])
    W     = float(sys.argv[2])
    H     = float(sys.argv[3])
    delta = float(sys.argv[4])
    gap   = float(sys.argv[5])
    L_DC  = float(sys.argv[6])

    r = coupling_coefficient(wl, W, H, delta, gap)
    if r is None:
        print("No guided mode."); sys.exit(1)
    kappa = r["kappa"]
    L_c = math.pi / (2 * kappa) if kappa > 0 else float("inf")
    L_50 = math.pi / (4 * kappa) if kappa > 0 else float("inf")
    P_cross   = math.sin(kappa * L_DC) ** 2
    P_through = math.cos(kappa * L_DC) ** 2
    L_cross_db = -10 * math.log10(P_cross)   if P_cross   > 0 else float("inf")
    L_through_db = -10 * math.log10(P_through) if P_through > 0 else float("inf")

    print(f"wavelength       = {wl:.4f} um")
    print(f"WG W x H         = {W:.4f} x {H:.4f} um")
    print(f"index contrast   = {delta:.4f} %")
    print(f"gap              = {gap:.4f} um")
    print(f"L_DC             = {L_DC:.4f} um   ({L_DC/1000:.4f} mm)")
    print(f"n_core           = {r['n_core']:.6f}")
    print(f"n_clad           = {r['n_clad']:.6f}")
    print(f"n_eff (total)    = {r['n_eff']:.6f}")
    print(f"kappa_x          = {r['kappa_x']:.6f} /um")
    print(f"gamma_x          = {r['gamma_x']:.6f} /um")
    print(f"beta             = {r['beta']:.6f} /um")
    print(f"kappa (coupling) = {kappa:.6e} /um  ({kappa*1000:.6f} /mm)")
    print(f"L_c (full)       = {L_c:.4f} um   ({L_c/1000:.4f} mm)")
    print(f"L_50 (3 dB)      = {L_50:.4f} um   ({L_50/1000:.4f} mm)")
    print(f"P_cross   @ L_DC = {P_cross:.6f}   ({L_cross_db:.3f} dB)")
    print(f"P_through @ L_DC = {P_through:.6f}   ({L_through_db:.3f} dB)")


if __name__ == "__main__":
    main()
