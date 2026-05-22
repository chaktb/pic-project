"""Bandpass filter from an asymmetric Mach-Zehnder interferometer.

The filter is a single asymmetric MZI: two identical directional couplers
(gap, L_DC) separated by two arms with a path-length difference Delta_L.
The cross port forms a bandpass centered where the arm phase Delta_phi = 0.

Design from a filter spec (waveguide W/H/Delta + DC gap given):
  * FSR   -> Delta_L = lambda_c^2 / (n_eff * FSR), snapped to an integer number
             of wavelengths so a cross-port peak sits exactly at lambda_c.
  * lambda_c (center) -> peak location and the wavelength axis.
  * passband -> coupler angle theta = kappa(lambda_c) * L_DC.
             For one MZI the cross response is
                 P_cross = sin^2(2 theta) * cos^2(Delta_phi / 2),
             so the absolute -3 dB (P > 0.5) bandwidth is tunable from 0 up to
             FSR/2 by under-coupling; the peak transmission sin^2(2 theta) is
             the price paid for a passband narrower than FSR/2.

Transfer function (matches mzi.py / MZI.java):
    P_bar = a^2 + b^2 - 2ab cos(Delta_phi),  a = cos(kL1)cos(kL2),
                                             b = sin(kL1)sin(kL2)
    P_cross = 1 - P_bar,   Delta_phi = (2 pi n_eff / lambda) * Delta_L

Usage:
    python band_pass_filter.py <W> <H> <delta_pct> <gap> <lam_c_um> <FSR_nm> <PB_nm> [span_xFSR=3] [N=600]
Example (1550 nm center, 1.6 nm FSR, 0.6 nm passband):
    python band_pass_filter.py 5.0 5.0 0.75 5.0 1.55 1.6 0.6
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


def kappa_neff(wl, W, H, delta, gap):
    """Returns (kappa per um, n_eff) at this wavelength, or None."""
    n_core = n_silica(wl)
    n_clad = n_core * (1 - delta / 100)
    v = slab_te_fundamental(wl, H, n_core, n_clad)
    if v is None:
        return None
    h = slab_te_fundamental(wl, W, v[2], n_clad)
    if h is None:
        return None
    kx, gx, n_eff = h
    k0 = 2 * math.pi / wl
    beta = n_eff * k0
    a = W / 2
    h2 = kx * kx
    kappa = 2 * h2 * gx * math.exp(-gx * gap) / (beta * (h2 + gx*gx) * (1 + gx * a))
    return kappa, n_eff


def mzi_power(wl, W, H, delta, gap1, L1, gap2, L2, dL):
    """Bar-port power of the 2-coupler MZI, or None."""
    kn1 = kappa_neff(wl, W, H, delta, gap1)
    kn2 = kappa_neff(wl, W, H, delta, gap2)
    if kn1 is None or kn2 is None:
        return None
    k1, n_eff = kn1
    k2, _ = kn2
    kL1 = k1 * L1
    kL2 = k2 * L2
    dphi = (2 * math.pi * n_eff / wl) * dL
    a_ = math.cos(kL1) * math.cos(kL2)
    b_ = math.sin(kL1) * math.sin(kL2)
    return a_*a_ + b_*b_ - 2*a_*b_*math.cos(dphi)


def design(W, H, delta, gap, lam_c, fsr_nm, pb_nm):
    """Solve for (Delta_L, L_DC) from the filter spec. Returns a dict or None."""
    kn = kappa_neff(lam_c, W, H, delta, gap)
    if kn is None:
        return None
    kappa_c, ne_c = kn
    fsr_um = fsr_nm / 1000.0
    pb_um = pb_nm / 1000.0

    # Delta_L from FSR, snapped to an integer number of wavelengths so that a
    # cross-port transmission peak (Delta_phi = 2*pi*m) sits exactly at lam_c.
    dL_target = lam_c * lam_c / (ne_c * fsr_um)
    m = max(1, round(ne_c * dL_target / lam_c))
    dL = m * lam_c / ne_c
    fsr_ach_um = lam_c * lam_c / (ne_c * dL)

    # Coupler angle from the requested absolute -3 dB passband (<= FSR/2).
    r = min(max(pb_um / fsr_ach_um, 1e-6), 0.5)
    theta = 0.5 * math.asin(min(1.0, math.sqrt(0.5) / math.cos(math.pi * r / 2)))
    l_dc = theta / kappa_c
    peak_t = math.sin(2 * theta) ** 2

    return {
        "kappa_c": kappa_c, "ne_c": ne_c, "m": m, "dL": dL,
        "fsr_ach_nm": fsr_ach_um * 1000.0, "theta": theta, "l_dc": l_dc,
        "peak_t": peak_t, "r": r, "clamped": pb_um / fsr_ach_um > 0.5,
    }


def main():
    if len(sys.argv) < 8:
        print("Usage: python band_pass_filter.py <W> <H> <delta_pct> <gap>"
              " <lam_c_um> <FSR_nm> <PB_nm> [span_xFSR=3] [N=600]")
        print("Example: python band_pass_filter.py 5.0 5.0 0.75 5.0 1.55 1.6 0.6")
        sys.exit(1)
    W = float(sys.argv[1])
    H = float(sys.argv[2])
    dl = float(sys.argv[3])
    gap = float(sys.argv[4])
    lam_c = float(sys.argv[5])
    fsr_nm = float(sys.argv[6])
    pb_nm = float(sys.argv[7])
    span = float(sys.argv[8]) if len(sys.argv) >= 9 else 3.0
    N = int(sys.argv[9]) if len(sys.argv) >= 10 else 600

    d = design(W, H, dl, gap, lam_c, fsr_nm, pb_nm)
    if d is None:
        print("# design failed: no guided mode for this waveguide/wavelength")
        sys.exit(1)

    l_dc = d["l_dc"]
    dL = d["dL"]
    fsr_ach_um = d["fsr_ach_nm"] / 1000.0

    print(f"# Bandpass (asymmetric MZI)  W={W} H={H} delta={dl}% gap={gap}")
    print(f"# spec: lambda_c={lam_c} um  FSR={fsr_nm} nm  passband={pb_nm} nm")
    print(f"#   n_eff(lam_c)   = {d['ne_c']:.6f}")
    print(f"#   kappa(lam_c)   = {d['kappa_c']:.6e} /um")
    print(f"#   Delta_L        = {dL:.4f} um   (m = {d['m']} wavelengths)")
    print(f"#   FSR achieved   = {d['fsr_ach_nm']:.4f} nm")
    print(f"#   coupler angle  = {d['theta']:.4f} rad ({math.degrees(d['theta']):.2f} deg)")
    print(f"#   L_DC (x2)      = {l_dc:.4f} um   gap = {gap} um")
    print(f"#   peak T (cross) = {d['peak_t']:.4f}")
    if d["clamped"]:
        print("#   NOTE: requested passband > FSR/2; clamped to FSR/2 (3 dB couplers)."
              " Use a cascade/lattice for narrower flat-top passbands.")
    print("# wavelength_um, P_bar, P_cross")

    wl_a = lam_c - 0.5 * span * fsr_ach_um
    wl_b = lam_c + 0.5 * span * fsr_ach_um
    for i in range(N):
        wl = wl_a + (wl_b - wl_a) * i / (N - 1)
        pb = mzi_power(wl, W, H, dl, gap, l_dc, gap, l_dc, dL)
        if pb is None:
            print(f"{wl:.6f}, NaN, NaN")
        else:
            print(f"{wl:.6f}, {pb:.6f}, {max(0.0, 1.0 - pb):.6f}")


if __name__ == "__main__":
    main()
