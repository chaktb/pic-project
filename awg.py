"""Arrayed Waveguide Grating (AWG) - N x N wavelength demultiplexer.

Two slab star couplers connected by an array of N_a waveguides with a constant
geometric path increment Delta_L. The grating order m fixes Delta_L = m*lam_c/n_eff_a;
the focal length f_l of the Rowland-circle slab sets the linear dispersion on the
output face so adjacent channels are separated by Delta_lam_ch.

Design equations (M.K. Smit, 1996; Takahashi, 1995):
    Delta_L = m * lam_c / n_eff_a                       (path increment)
    FSR_lam = lam_c^2 / (n_g * Delta_L)                 (free spectral range)
    f_l     = (n_eff_s * d_aw * d_io * n_eff_a) / (m * n_g * Delta_lam_ch)
    theta_i = (i - (N_ch+1)/2) * d_io / f_l             (port angle)

Output transmission at port i:
    T_i(lam) = | sum_k A_k * exp(j*(k-c)*phi_i(lam)) |^2 / (sum_k A_k)^2
    phi_i(lam) = beta(lam)*Delta_L + 2*pi*n_eff_s*d_aw*theta_i / lam
    beta = 2*pi*n_eff_a(lam) / lam

A_k is the coupling profile at the slab/array junction: uniform or Gaussian
(width parameter w_g typ. 0.3-0.4 of N_a).

Refractive indices follow Sellmeier silica (Malitson 1965) + Effective Index
Method (TE-TE) for the channel waveguide; n_eff_s is the vertical slab mode
index in the star coupler.

Usage:
    python awg.py <W> <H> <delta_pct> <lam_c_um> <dlam_ch_nm> <N_ch> <N_a> <d_aw_um> <d_io_um> [m=auto] [profile=gauss|uniform] [w_g=0.35] [span=1.0] [N=900]
Example (1550 nm grid, 0.8 nm spacing, 1x4 AWG):
    python awg.py 6.0 6.0 0.75 1.55 0.8 4 30 15.0 20.0
"""

import math
import sys


def n_silica(wl):
    B1, B2, B3 = 0.6961663, 0.4079426, 0.8974794
    C1, C2, C3 = 0.0684043**2, 0.1162414**2, 9.896161**2
    l2 = wl * wl
    return math.sqrt(1 + B1*l2/(l2-C1) + B2*l2/(l2-C2) + B3*l2/(l2-C3))


def slab_te_fundamental(wl, d, n_core, n_clad):
    """Fundamental TE mode of a symmetric slab (kappa, n_eff). Returns None if cut off."""
    if n_core <= n_clad:
        return None
    k0 = 2*math.pi / wl
    a = d / 2.0
    V = k0 * a * math.sqrt(n_core*n_core - n_clad*n_clad)
    upper = min(math.pi/2 - 1e-9, V - 1e-12)
    if upper <= 0:
        return None
    def f(u):
        w = math.sqrt(max(V*V - u*u, 0.0))
        return u * math.tan(u) - w
    lo, hi = 1e-9, upper
    flo = f(lo); fhi = f(hi)
    if flo * fhi > 0:
        return None
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        fmid = f(mid)
        if flo * fmid <= 0:
            hi = mid; fhi = fmid
        else:
            lo = mid; flo = fmid
        if hi - lo < 1e-13:
            break
    u = 0.5 * (lo + hi)
    return u/a, math.sqrt(n_core*n_core - (u/a/k0)**2)


def eim(wl, W, H, delta_pct):
    """Effective Index Method (TE-TE) for the channel waveguide.

    Returns (n_eff_a, n_eff_s):
        n_eff_a = channel mode index (used as beta source for the array waveguide)
        n_eff_s = vertical slab mode index (used for the star coupler region)
    """
    n_clad = n_silica(wl)
    n_core = n_clad / (1 - delta_pct/100.0)
    v = slab_te_fundamental(wl, H, n_core, n_clad)
    if v is None:
        return None
    _, n_eff_s = v
    h = slab_te_fundamental(wl, W, n_eff_s, n_clad)
    if h is None:
        return None
    _, n_eff_a = h
    return n_eff_a, n_eff_s


def group_index(wl, W, H, delta_pct):
    """n_g = n_eff_a - lam * d(n_eff_a)/d(lam) via central difference."""
    dwl = wl * 1e-4
    a = eim(wl - dwl, W, H, delta_pct)
    b = eim(wl + dwl, W, H, delta_pct)
    c = eim(wl, W, H, delta_pct)
    if a is None or b is None or c is None:
        return None
    dn_dwl = (b[0] - a[0]) / (2*dwl)
    return c[0] - wl * dn_dwl


def design(W, H, delta_pct, lam_c, dlam_ch_um, n_ch, n_a, d_aw, d_io, m_pref=0):
    """Return the AWG design parameters as a dict (all lengths in um)."""
    n_eff_a, n_eff_s = eim(lam_c, W, H, delta_pct)
    n_g = group_index(lam_c, W, H, delta_pct) or n_eff_a

    # If user did not specify the phase order, pick m such that FSR >= (N_ch+2)*Dlam_ch.
    # m = round(lam_c * n_eff_a / (n_g * FSR_target))
    if m_pref <= 0:
        fsr_target = (n_ch + 2) * dlam_ch_um
        m = max(1, round(lam_c * n_eff_a / (n_g * fsr_target)))
    else:
        m = int(m_pref)

    dL = m * lam_c / n_eff_a
    fsr = lam_c * lam_c / (n_g * dL)
    fl = (n_eff_s * d_aw * d_io * n_eff_a) / (m * n_g * dlam_ch_um)
    sweep = (n_ch - 1) * d_io / fl
    return dict(n_eff_a=n_eff_a, n_eff_s=n_eff_s, n_g=n_g,
                m=m, dL=dL, fsr=fsr, fl=fl, sweep=sweep)


def coupling_amp(n_a, profile="gauss", w_g=0.35):
    A = [0.0] * (n_a + 1)
    c = 0.5 * (n_a + 1)
    for k in range(1, n_a + 1):
        if profile == "gauss":
            x = (k - c) / (w_g * n_a)
            A[k] = math.exp(-x*x)
        else:
            A[k] = 1.0
    return A


def transmission(wl, port, p, d, A):
    """Power transmission at output `port` (1..N_ch) at wavelength wl [um]."""
    res = eim(wl, p["W"], p["H"], p["delta_pct"])
    if res is None:
        return float("nan")
    n_eff_a, n_eff_s = res
    beta = 2*math.pi * n_eff_a / wl
    cch = 0.5 * (p["n_ch"] + 1)
    theta = (port - cch) * p["d_io"] / d["fl"]
    dphi = beta * d["dL"] + 2*math.pi * n_eff_s * p["d_aw"] * theta / wl
    c = 0.5 * (p["n_a"] + 1)
    re = im = a_sum = 0.0
    for k in range(1, p["n_a"] + 1):
        phi = (k - c) * dphi
        re += A[k] * math.cos(phi)
        im += A[k] * math.sin(phi)
        a_sum += A[k]
    return (re*re + im*im) / (a_sum * a_sum)


def main(argv):
    if len(argv) < 10:
        print(__doc__)
        sys.exit(1)
    W = float(argv[1]); H = float(argv[2]); delta_pct = float(argv[3])
    lam_c = float(argv[4]); dlam_ch_nm = float(argv[5])
    n_ch = int(argv[6]); n_a = int(argv[7])
    d_aw = float(argv[8]); d_io = float(argv[9])
    m_pref = int(argv[10]) if len(argv) > 10 else 0
    profile = argv[11] if len(argv) > 11 else "gauss"
    w_g = float(argv[12]) if len(argv) > 12 else 0.35
    span = float(argv[13]) if len(argv) > 13 else 1.0
    N = int(argv[14]) if len(argv) > 14 else 900

    dlam_ch_um = dlam_ch_nm * 1e-3
    p = dict(W=W, H=H, delta_pct=delta_pct, lam_c=lam_c, d_aw=d_aw, d_io=d_io,
             n_ch=n_ch, n_a=n_a)
    d = design(W, H, delta_pct, lam_c, dlam_ch_um, n_ch, n_a, d_aw, d_io, m_pref)
    A = coupling_amp(n_a, profile, w_g)

    print(f"n_eff_a = {d['n_eff_a']:.6f}")
    print(f"n_eff_s = {d['n_eff_s']:.6f}")
    print(f"n_g     = {d['n_g']:.6f}")
    print(f"m       = {d['m']}")
    print(f"Delta_L = {d['dL']:.4f} um")
    print(f"FSR     = {d['fsr']*1000:.4f} nm")
    print(f"f_l     = {d['fl']*1e-3:.4f} mm")
    print(f"sweep   = {math.degrees(d['sweep']):.3f} deg")

    print()
    print("# wavelength_um  " + "  ".join(f"port{i}_dB" for i in range(1, n_ch+1)))
    half = 0.5 * span * d["fsr"]
    for n in range(N):
        wl = lam_c - half + 2*half * n / (N - 1)
        row = [f"{wl:.6f}"]
        for i in range(1, n_ch+1):
            T = transmission(wl, i, p, d, A)
            dB = 10*math.log10(T) if T > 1e-12 else -120.0
            row.append(f"{dB:8.3f}")
        print("  ".join(row))


if __name__ == "__main__":
    main(sys.argv)
