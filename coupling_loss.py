"""Butt-coupling loss between SMF and a rectangular channel waveguide.

Method:
  - SMF fundamental mode: circular Gaussian E_SMF(x,y) = exp(-(x^2+y^2)/w_smf^2),
    w_smf = MFD/2.
  - Channel waveguide: separable mode E_WG(x,y) = X(x)*Y(y) from EIM (TE-TE).
    Each 1D profile is the symmetric-slab fundamental: cos(kx*r) inside the core,
    cos(kx*a)*exp(-gamma*(|r|-a)) outside.
  - Overlap factorizes:  eta = eta_x * eta_y, with
    eta_dir = (integral G * F)^2 / (integral G^2 * integral F^2).
  - Coupling loss (dB) = -10 log10(eta).

Sellmeier (Malitson 1965) gives n_core for fused silica.
"""

import math
import sys


def n_silica(wl):
    B1, B2, B3 = 0.6961663, 0.4079426, 0.8974794
    C1, C2, C3 = 0.0684043**2, 0.1162414**2, 9.896161**2
    l2 = wl * wl
    n2 = 1 + B1*l2/(l2-C1) + B2*l2/(l2-C2) + B3*l2/(l2-C3)
    return math.sqrt(n2)


def slab_te_fundamental(wl, d, n_core, n_clad):
    """Returns (kappa, gamma, n_eff) for the symmetric-slab fundamental TE mode.

    kappa: transverse wavenumber inside the core
    gamma: decay constant in the cladding
    n_eff: modal effective index
    Returns None if no guided mode.
    """
    if n_core <= n_clad:
        return None
    k0 = 2 * math.pi / wl
    a = d / 2.0
    V = k0 * a * math.sqrt(n_core**2 - n_clad**2)
    if V <= 0:
        return None
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
    kappa = u / a
    gamma = w_v / a
    n_eff = math.sqrt(n_core**2 - (kappa / k0)**2)
    return kappa, gamma, n_eff


def slab_field(r, kappa, gamma, a):
    if abs(r) <= a:
        return math.cos(kappa * r)
    return math.cos(kappa * a) * math.exp(-gamma * (abs(r) - a))


def gaussian(r, w):
    return math.exp(-(r*r) / (w*w))


def integrate_trap(f, lo, hi, n):
    h = (hi - lo) / n
    s = 0.5 * (f(lo) + f(hi))
    for i in range(1, n):
        s += f(lo + i * h)
    return s * h


def coupling_loss(wl, mfd_smf, W, H, n_core, n_clad):
    vmode = slab_te_fundamental(wl, H, n_core, n_clad)
    if vmode is None:
        return None
    kappa_y, gamma_y, n_eff_v = vmode
    hmode = slab_te_fundamental(wl, W, n_eff_v, n_clad)
    if hmode is None:
        return None
    kappa_x, gamma_x, _ = hmode
    a_x, a_y = W / 2, H / 2
    w_smf = mfd_smf / 2

    L = max(6 * w_smf, 6 * W, 6 * H)
    N = 4000

    def overlap(field, half):
        num = integrate_trap(lambda r: gaussian(r, w_smf) * field(r), -L, L, N)
        den_g = integrate_trap(lambda r: gaussian(r, w_smf) ** 2, -L, L, N)
        den_f = integrate_trap(lambda r: field(r) ** 2, -L, L, N)
        if den_g <= 0 or den_f <= 0:
            return 0.0
        return (num * num) / (den_g * den_f)

    eta_x = overlap(lambda x: slab_field(x, kappa_x, gamma_x, a_x), a_x)
    eta_y = overlap(lambda y: slab_field(y, kappa_y, gamma_y, a_y), a_y)
    eta = eta_x * eta_y
    loss_db = -10 * math.log10(eta) if eta > 0 else float("inf")
    return {
        "n_eff_vertical": n_eff_v,
        "kappa_x": kappa_x, "gamma_x": gamma_x,
        "kappa_y": kappa_y, "gamma_y": gamma_y,
        "eta_x": eta_x, "eta_y": eta_y,
        "eta": eta, "loss_db": loss_db,
    }


def main():
    if len(sys.argv) != 6:
        print("Usage: python coupling_loss.py <wavelength_um> <MFD_smf_um> <W_um> <H_um> <delta_percent>")
        print("Example: python coupling_loss.py 1.55 10.4 5.0 5.0 0.75")
        sys.exit(1)
    wl, mfd, W, H, delta = (float(x) for x in sys.argv[1:6])
    n_core = n_silica(wl)
    n_clad = n_core * (1 - delta / 100)

    r = coupling_loss(wl, mfd, W, H, n_core, n_clad)
    if r is None:
        print("No guided mode under given parameters.")
        sys.exit(1)

    print(f"wavelength            = {wl:.4f} um")
    print(f"SMF MFD               = {mfd:.4f} um  (w_smf = {mfd/2:.4f} um)")
    print(f"WG core W x H         = {W:.4f} x {H:.4f} um")
    print(f"index contrast        = {delta:.4f} %")
    print(f"n_core (Sellmeier)    = {n_core:.6f}")
    print(f"n_clad                = {n_clad:.6f}")
    print(f"n_eff (vertical slab) = {r['n_eff_vertical']:.6f}")
    print(f"eta_x                 = {r['eta_x']:.6f}")
    print(f"eta_y                 = {r['eta_y']:.6f}")
    print(f"eta (total)           = {r['eta']:.6f}")
    print(f"coupling loss         = {r['loss_db']:.3f} dB")


if __name__ == "__main__":
    main()
