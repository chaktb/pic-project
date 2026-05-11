"""Mach-Zehnder Interferometer (MZI) spectral response.

Two directional couplers (each: gap, L_DC) with two arms whose path lengths
differ by Delta_L. The transfer function:

    E_bar  (out port 1) = cos(k1)cos(k2) - sin(k1)sin(k2) * exp(-j*dphi)
    E_cross(out port 2) = -j [ sin(k2)cos(k1) + cos(k2)sin(k1) * exp(-j*dphi) ]

with k_i = kappa(lambda) * L_DC_i and dphi = beta(lambda) * Delta_L.
For two identical 3 dB couplers (k=pi/4):
    P_bar   = sin^2(dphi/2)
    P_cross = cos^2(dphi/2)

FSR ~ lambda^2 / (n_eff * Delta_L) — the channel spacing.
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


def kappa_and_neff(wl, W, H, delta, gap):
    """Returns (kappa per um, n_eff) at this wavelength, or (None, None)."""
    n_core = n_silica(wl)
    n_clad = n_core * (1 - delta / 100)
    vmode = slab_te_fundamental(wl, H, n_core, n_clad)
    if vmode is None:
        return None, None
    _, _, n_eff_v = vmode
    hmode = slab_te_fundamental(wl, W, n_eff_v, n_clad)
    if hmode is None:
        return None, None
    kappa_x, gamma_x, n_eff = hmode
    k0 = 2 * math.pi / wl
    beta = n_eff * k0
    a = W / 2
    h2 = kappa_x * kappa_x
    kappa = 2 * h2 * gamma_x * math.exp(-gamma_x * gap) / (beta * (h2 + gamma_x*gamma_x) * (1 + gamma_x * a))
    return kappa, n_eff


def mzi_outputs(wl, W, H, delta, gap, L_DC, dL):
    kappa, n_eff = kappa_and_neff(wl, W, H, delta, gap)
    if kappa is None:
        return None
    k1 = k2 = kappa * L_DC
    dphi = (2 * math.pi * n_eff / wl) * dL
    a_ = math.cos(k1) * math.cos(k2)
    b_ = math.sin(k1) * math.sin(k2)
    P_bar = a_*a_ + b_*b_ - 2*a_*b_*math.cos(dphi)
    P_cross = max(0.0, 1.0 - P_bar)
    return {"kappa": kappa, "n_eff": n_eff, "kL": k1, "dphi": dphi,
            "P_bar": P_bar, "P_cross": P_cross}


def main():
    if len(sys.argv) not in (10, 11):
        print("Usage: python mzi.py <wl_start> <wl_end> <W> <H> <delta_pct> <gap> <L_DC> <Delta_L> [N=400]")
        print("Example: python mzi.py 1.50 1.60 5.0 5.0 0.75 5.0 300 100")
        sys.exit(1)
    wl_a = float(sys.argv[1])
    wl_b = float(sys.argv[2])
    W    = float(sys.argv[3])
    H    = float(sys.argv[4])
    dl   = float(sys.argv[5])
    gap  = float(sys.argv[6])
    L_DC = float(sys.argv[7])
    dL   = float(sys.argv[8])
    N    = int(sys.argv[9]) if len(sys.argv) >= 10 else 400

    wlc = 0.5 * (wl_a + wl_b)
    rc = mzi_outputs(wlc, W, H, dl, gap, L_DC, dL)
    if rc:
        fsr_um = wlc * wlc / (rc["n_eff"] * max(dL, 1e-12))
        print(f"# MZI W={W} H={H} delta={dl}% gap={gap} L_DC={L_DC} dL={dL}")
        print(f"# at lambda_center={wlc:.4f}: kappa={rc['kappa']:.4e}/um  kL={rc['kL']:.4f}  n_eff={rc['n_eff']:.6f}")
        print(f"# approximate FSR ~ {fsr_um*1000:.4f} nm")
    print("# wavelength_um, P_bar, P_cross")
    for i in range(N):
        wl = wl_a + (wl_b - wl_a) * i / (N - 1)
        r = mzi_outputs(wl, W, H, dl, gap, L_DC, dL)
        if r is None:
            print(f"{wl:.6f}, NaN, NaN")
        else:
            print(f"{wl:.6f}, {r['P_bar']:.6f}, {r['P_cross']:.6f}")


if __name__ == "__main__":
    main()
