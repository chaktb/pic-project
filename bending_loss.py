"""Bending (radiation) loss for a rectangular silica channel waveguide.

Assumes in-plane bending (radius R measured in the W-direction). Vertical mode
solved first for n_eff_v; horizontal slab (W) then gives kappa_x, gamma_x,
beta = k0 * n_eff_total.

Per-unit-length radiation loss (Marcuse-style approximation):

    alpha(R) = (gamma/2) * exp(-2*gamma*a) * exp(-(2*gamma^3 / (3*beta^2)) * R)

Total bend loss for an arc of angle theta (radians):

    Loss_dB(R, theta) = 4.343 * alpha(R) * R * theta

Sellmeier (Malitson 1965) gives n_core for fused silica.
"""

import math
import sys


def n_silica(wl):
    B1, B2, B3 = 0.6961663, 0.4079426, 0.8974794
    C1, C2, C3 = 0.0684043**2, 0.1162414**2, 9.896161**2
    l2 = wl * wl
    return math.sqrt(1 + B1*l2/(l2-C1) + B2*l2/(l2-C2) + B3*l2/(l2-C3))


def slab_te_fundamental(wl, d, n_core, n_clad):
    """Returns (kappa, gamma, n_eff) or None."""
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


def bending_alpha_per_um(R_um, gamma_um, beta_um, a_um):
    """Per-micron radiation loss coefficient at bend radius R (all in microns)."""
    A = (gamma_um / 2) * math.exp(-2 * gamma_um * a_um)
    B = (2 * gamma_um**3) / (3 * beta_um**2)
    return A * math.exp(-B * R_um)


def bending_loss_db(R_um, theta_rad, gamma_um, beta_um, a_um):
    alpha = bending_alpha_per_um(R_um, gamma_um, beta_um, a_um)
    return 4.342944819 * alpha * R_um * theta_rad


def critical_R_for_loss(target_db, theta_rad, gamma_um, beta_um, a_um):
    """Smallest R [um] in the LARGE-R monotonic region where bend loss <= target_db.

    Loss(R) = 4.343 * A * R * theta * exp(-B*R) is non-monotonic with peak at
    R_peak = 1/B. The design-relevant root is on the right side (R > R_peak):
    the smallest R for which loss is acceptable as R grows.
    """
    B = (2 * gamma_um**3) / (3 * beta_um**2)
    if B <= 0:
        return None
    R_peak = 1.0 / B
    # If the peak itself is below the target, any R >= 0 works -> return 0
    if bending_loss_db(R_peak, theta_rad, gamma_um, beta_um, a_um) <= target_db:
        return 0.0
    # Find R_hi above the peak with loss <= target
    R_hi = R_peak * 2
    for _ in range(80):
        if bending_loss_db(R_hi, theta_rad, gamma_um, beta_um, a_um) <= target_db:
            break
        R_hi *= 2
        if R_hi > 1e12:
            return None
    R_lo = R_peak  # loss above target on this side
    for _ in range(120):
        mid = 0.5 * (R_lo + R_hi)
        if bending_loss_db(mid, theta_rad, gamma_um, beta_um, a_um) <= target_db:
            R_hi = mid
        else:
            R_lo = mid
        if R_hi - R_lo < 1e-3:
            break
    return R_hi


def main():
    if len(sys.argv) not in (5, 6, 7):
        print("Usage: python bending_loss.py <wl_um> <W_um> <H_um> <delta_pct> [theta_deg=90] [R_mm=5]")
        print("Example: python bending_loss.py 1.55 5.0 5.0 0.75 90 5")
        sys.exit(1)
    wl    = float(sys.argv[1])
    W     = float(sys.argv[2])
    H     = float(sys.argv[3])
    delta = float(sys.argv[4])
    theta_deg = float(sys.argv[5]) if len(sys.argv) >= 6 else 90.0
    R_mm  = float(sys.argv[6]) if len(sys.argv) >= 7 else 5.0
    theta = math.radians(theta_deg)
    R_um = R_mm * 1000

    n_core = n_silica(wl)
    n_clad = n_core * (1 - delta / 100)
    vmode = slab_te_fundamental(wl, H, n_core, n_clad)
    if vmode is None:
        print("No vertical mode."); sys.exit(1)
    _, _, n_eff_v = vmode
    hmode = slab_te_fundamental(wl, W, n_eff_v, n_clad)
    if hmode is None:
        print("No horizontal mode."); sys.exit(1)
    kappa_x, gamma_x, n_eff = hmode
    k0 = 2 * math.pi / wl
    beta = n_eff * k0
    a_x = W / 2

    alpha_um = bending_alpha_per_um(R_um, gamma_x, beta, a_x)
    loss_db = bending_loss_db(R_um, theta, gamma_x, beta, a_x)
    Rc_01 = critical_R_for_loss(0.1, theta, gamma_x, beta, a_x)
    Rc_1  = critical_R_for_loss(1.0, theta, gamma_x, beta, a_x)

    print(f"wavelength       = {wl:.4f} um")
    print(f"W x H            = {W:.4f} x {H:.4f} um")
    print(f"index contrast   = {delta:.4f} %")
    print(f"bend angle       = {theta_deg:.2f} deg")
    print(f"R (point readout)= {R_mm:.4f} mm")
    print(f"n_core           = {n_core:.6f}")
    print(f"n_clad           = {n_clad:.6f}")
    print(f"n_eff (total)    = {n_eff:.6f}")
    print(f"gamma_x          = {gamma_x:.6f} /um")
    print(f"beta             = {beta:.6f} /um")
    print(f"alpha (at R)     = {alpha_um:.4e} /um  ({alpha_um*1e3:.4e} /mm)")
    print(f"loss @ R         = {loss_db:.4f} dB / {theta_deg:g}deg bend")
    if Rc_01 is not None:
        print(f"R_critical(0.1dB)= {Rc_01/1000:.4f} mm")
    if Rc_1 is not None:
        print(f"R_critical(1dB)  = {Rc_1/1000:.4f} mm")


if __name__ == "__main__":
    main()
