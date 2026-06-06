"""Thermally tunable bandpass filter (asymmetric MZI + NiCr thin-film heater).

The filter is an asymmetric MZI: two identical, *tunable* directional couplers
(each a balanced-MZI coupler under a NiCr heater) separated by arms with a
fixed path-length difference Delta_L. Because Delta_L is fixed and snapped to an
integer number of wavelengths at lambda_c, the cross-port transmission peak
stays locked at lambda_c (the center wavelength does NOT move). The NiCr heater
changes only the coupler angle theta, and hence only the passband WIDTH.

Heater chain (NiCr thin film over the couplers):
    R   = rho_NiCr * L_h / (w_h * t_h)            [film resistance]
    P   = V^2 / R                                  [dissipated power]
    dT  = R_th * P                                 [temperature rise]
    dn  = (dn/dT) * dT                             [thermo-optic index change]
    d_theta = pi * dn * L_h / lambda_c             [extra coupling angle]
    theta(T) = theta_base + d_theta

Cross-port response (matches band_pass_filter.py / mzi):
    P_cross = sin^2(2 theta) * cos^2(Delta_phi/2),  Delta_phi = (2 pi n_eff/lambda) Delta_L
The absolute -3 dB passband width is
    PB = FSR * (2/pi) * acos( sqrt(0.5 / sin^2(2 theta)) )   (0 if peak <= 0.5)
so heating widens the passband toward FSR/2 (at theta = pi/4) while lambda_c
stays fixed. Peak transmission sin^2(2 theta) is the single-stage trade-off.

Usage:
    python tunable_bandpass.py <W> <H> <delta_pct> <gap> <lam_c_um> <FSR_nm> <base_PB_nm> <V> [span_xFSR=3] [N=600]
Example (1550 nm fixed center, 1.6 nm FSR, 0.3 nm base passband, 5 V):
    python tunable_bandpass.py 5.0 5.0 0.75 5.0 1.55 1.6 0.3 5.0
"""

import math
import sys

# NiCr / silica defaults (overridable below if you edit the constants).
RHO_NICR = 1.1e-6     # Ohm*m  (~110 uOhm*cm)
DNDT = 1.0e-5         # 1/K    (fused silica thermo-optic coefficient)
R_TH = 2000.0         # K/W    (lumped heater-to-waveguide thermal resistance)
W_H = 4.0             # um     (heater width)
T_H = 0.1             # um     (heater film thickness, 100 nm)
# Heater length is taken equal to the coupler length L_DC (it sits over it).


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
    n_clad = n_silica(wl)
    n_core = n_clad / (1 - delta / 100)
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


def mzi_power(wl, W, H, delta, gap1, L1, gap2, L2, dL, dtheta=0.0):
    """Bar-port power of the 2-coupler MZI, or None.

    dtheta is the heater-induced extra coupling angle added to each (tunable)
    coupler, so the effective angle is kappa(wl)*L + dtheta.
    """
    kn1 = kappa_neff(wl, W, H, delta, gap1)
    kn2 = kappa_neff(wl, W, H, delta, gap2)
    if kn1 is None or kn2 is None:
        return None
    k1, n_eff = kn1
    k2, _ = kn2
    kL1 = k1 * L1 + dtheta
    kL2 = k2 * L2 + dtheta
    dphi = (2 * math.pi * n_eff / wl) * dL
    a_ = math.cos(kL1) * math.cos(kL2)
    b_ = math.sin(kL1) * math.sin(kL2)
    return a_*a_ + b_*b_ - 2*a_*b_*math.cos(dphi)


def pb_from_theta(theta, fsr):
    """Absolute -3 dB passband width (same units as fsr) for coupler angle theta."""
    peak = math.sin(2 * theta) ** 2
    if peak <= 0.5:
        return 0.0
    return fsr * (2.0 / math.pi) * math.acos(math.sqrt(0.5 / peak))


def heater(volt, l_h_um):
    """NiCr heater electricals/thermals. Returns (R, P, dT, dn)."""
    l_m = l_h_um * 1e-6
    w_m = W_H * 1e-6
    t_m = T_H * 1e-6
    R = RHO_NICR * l_m / (w_m * t_m)
    P = volt * volt / R if R > 0 else 0.0
    dT = R_TH * P
    dn = DNDT * dT
    return R, P, dT, dn


def design(W, H, delta, gap, lam_c, fsr_nm, base_pb_nm, volt):
    """Design + heater state. Returns a dict, or None if no guided mode."""
    kn = kappa_neff(lam_c, W, H, delta, gap)
    if kn is None:
        return None
    kappa_c, ne_c = kn
    fsr_um = fsr_nm / 1000.0
    base_pb_um = base_pb_nm / 1000.0

    # Fixed center: snap Delta_L to an integer number of wavelengths at lam_c.
    dL_target = lam_c * lam_c / (ne_c * fsr_um)
    m = max(1, round(ne_c * dL_target / lam_c))
    dL = m * lam_c / ne_c
    fsr_ach_um = lam_c * lam_c / (ne_c * dL)

    # Base coupler angle from the base (unpowered) passband (<= FSR/2).
    r0 = min(max(base_pb_um / fsr_ach_um, 1e-6), 0.5)
    theta_base = 0.5 * math.asin(min(1.0, math.sqrt(0.5) / math.cos(math.pi * r0 / 2)))
    l_dc = theta_base / kappa_c

    # Heater (length = coupler length) tunes theta, hence width only.
    R, P, dT, dn = heater(volt, l_dc)
    d_theta = math.pi * dn * l_dc / lam_c
    theta_t = theta_base + d_theta
    peak_t = math.sin(2 * theta_t) ** 2

    return {
        "kappa_c": kappa_c, "ne_c": ne_c, "m": m, "dL": dL,
        "fsr_ach_um": fsr_ach_um, "fsr_ach_nm": fsr_ach_um * 1000.0,
        "theta_base": theta_base, "l_dc": l_dc,
        "R": R, "P": P, "dT": dT, "dn": dn, "d_theta": d_theta,
        "theta_t": theta_t, "peak_t": peak_t,
        "pb_base_nm": pb_from_theta(theta_base, fsr_ach_um) * 1000.0,
        "pb_tuned_nm": pb_from_theta(theta_t, fsr_ach_um) * 1000.0,
    }


def main():
    if len(sys.argv) < 9:
        print("Usage: python tunable_bandpass.py <W> <H> <delta_pct> <gap>"
              " <lam_c_um> <FSR_nm> <base_PB_nm> <V> [span_xFSR=3] [N=600]")
        print("Example: python tunable_bandpass.py 5.0 5.0 0.75 5.0 1.55 1.6 0.3 5.0")
        sys.exit(1)
    W = float(sys.argv[1])
    H = float(sys.argv[2])
    dl = float(sys.argv[3])
    gap = float(sys.argv[4])
    lam_c = float(sys.argv[5])
    fsr_nm = float(sys.argv[6])
    base_pb_nm = float(sys.argv[7])
    volt = float(sys.argv[8])
    span = float(sys.argv[9]) if len(sys.argv) >= 10 else 3.0
    N = int(sys.argv[10]) if len(sys.argv) >= 11 else 600

    d = design(W, H, dl, gap, lam_c, fsr_nm, base_pb_nm, volt)
    if d is None:
        print("# design failed: no guided mode for this waveguide/wavelength")
        sys.exit(1)

    l_dc = d["l_dc"]
    dL = d["dL"]
    fsr_ach_um = d["fsr_ach_um"]

    print(f"# Tunable bandpass (asymmetric MZI + NiCr heater)  W={W} H={H} delta={dl}% gap={gap}")
    print(f"# spec: lambda_c={lam_c} um (FIXED)  FSR={fsr_nm} nm  base_PB={base_pb_nm} nm  V={volt} V")
    print(f"#   n_eff(lam_c)   = {d['ne_c']:.6f}")
    print(f"#   Delta_L        = {dL:.4f} um   (m = {d['m']} wavelengths, center locked)")
    print(f"#   FSR achieved   = {d['fsr_ach_nm']:.4f} nm")
    print(f"#   L_DC (x2)      = {l_dc:.4f} um   (heater length = L_DC)")
    print(f"#   theta_base     = {math.degrees(d['theta_base']):.3f} deg   base PB = {d['pb_base_nm']:.4f} nm")
    print("#   -- NiCr heater --")
    print(f"#   R              = {d['R']:.2f} Ohm")
    print(f"#   P = V^2/R      = {d['P']*1e3:.4f} mW")
    print(f"#   dT             = {d['dT']:.3f} K")
    print(f"#   dn (T-O)       = {d['dn']:.3e}")
    print(f"#   d_theta        = {math.degrees(d['d_theta']):.3f} deg")
    print(f"#   theta(T)       = {math.degrees(d['theta_t']):.3f} deg")
    print(f"#   peak T (cross) = {d['peak_t']:.4f}")
    print(f"#   PB tuned       = {d['pb_tuned_nm']:.4f} nm   (center unchanged at {lam_c} um)")
    print("# wavelength_um, P_bar, P_cross")

    wl_a = lam_c - 0.5 * span * fsr_ach_um
    wl_b = lam_c + 0.5 * span * fsr_ach_um
    for i in range(N):
        wl = wl_a + (wl_b - wl_a) * i / (N - 1)
        pb = mzi_power(wl, W, H, dl, gap, l_dc, gap, l_dc, dL, d["d_theta"])
        if pb is None:
            print(f"{wl:.6f}, NaN, NaN")
        else:
            print(f"{wl:.6f}, {pb:.6f}, {max(0.0, 1.0 - pb):.6f}")


if __name__ == "__main__":
    main()
