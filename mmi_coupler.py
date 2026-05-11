"""M x N MMI coupler dimension design.

Approach (general interference, Soldano & Pennings 1995):
  - pitch  = W_wg + gap                       (port pitch at MMI edges)
  - W_MMI  = max(M, N) * pitch                (MMI slab width)
  - W_e    = W_MMI + (lambda/pi) * (n_clad/n_core)^(2*sigma) / sqrt(n_core^2 - n_clad^2)
            (Goos-Hanchen penetration; sigma=0 for TE, 1 for TM)
  - L_pi   = 4 * n_core * W_e^2 / (3 * lambda)
  - L_MMI  = 3 * L_pi / N                     (general interference, N images at output)

Sellmeier (Malitson 1965) for n_core(lambda).
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


def mmi_design(wl, W_wg, H, delta, M, N, gap, sigma=0):
    n_core = n_silica(wl)
    n_clad = n_core * (1 - delta / 100.0)
    vmode = slab_te_fundamental(wl, H, n_core, n_clad)
    n_eff = vmode[2] if vmode else n_core

    Nmax = max(M, N)
    pitch = W_wg + gap
    W_MMI = Nmax * pitch
    W_e = W_MMI + (wl / math.pi) * ((n_clad / n_core) ** (2 * sigma)) / math.sqrt(n_core**2 - n_clad**2)
    L_pi = 4 * n_core * W_e * W_e / (3 * wl)
    L_MMI = 3 * L_pi / N

    # Port positions (centered at y=0). Same pitch on input and output sides.
    y0 = -(Nmax - 1) / 2.0 * pitch
    in_ports  = [y0 + (Nmax - M) / 2.0 * pitch + i * pitch for i in range(M)]
    out_ports = [y0 + (Nmax - N) / 2.0 * pitch + i * pitch for i in range(N)]

    return {
        "n_core": n_core, "n_clad": n_clad, "n_eff": n_eff,
        "pitch": pitch, "W_MMI": W_MMI, "W_e": W_e,
        "L_pi": L_pi, "L_MMI": L_MMI,
        "in_ports": in_ports, "out_ports": out_ports,
    }


def main():
    if len(sys.argv) != 8:
        print("Usage: python mmi_coupler.py <wl_um> <W_wg_um> <H_um> <delta_pct> <M> <N> <gap_um>")
        print("Example: python mmi_coupler.py 1.55 5.0 5.0 0.75 2 2 5.0")
        sys.exit(1)
    wl    = float(sys.argv[1])
    W_wg  = float(sys.argv[2])
    H     = float(sys.argv[3])
    delta = float(sys.argv[4])
    M     = int(sys.argv[5])
    N     = int(sys.argv[6])
    gap   = float(sys.argv[7])

    r = mmi_design(wl, W_wg, H, delta, M, N, gap)
    print(f"wavelength      = {wl:.4f} um")
    print(f"WG W x H        = {W_wg:.4f} x {H:.4f} um")
    print(f"index contrast  = {delta:.4f} %")
    print(f"M x N           = {M} x {N}")
    print(f"gap             = {gap:.4f} um")
    print(f"n_core          = {r['n_core']:.6f}")
    print(f"n_clad          = {r['n_clad']:.6f}")
    print(f"n_eff (vertical)= {r['n_eff']:.6f}")
    print(f"pitch           = {r['pitch']:.4f} um")
    print(f"W_MMI           = {r['W_MMI']:.4f} um")
    print(f"W_e (eff width) = {r['W_e']:.4f} um")
    print(f"L_pi (beat)     = {r['L_pi']:.4f} um   ({r['L_pi']/1000:.4f} mm)")
    print(f"L_MMI           = {r['L_MMI']:.4f} um   ({r['L_MMI']/1000:.4f} mm)")
    print(f"input  port y   = {[round(y,3) for y in r['in_ports']]}")
    print(f"output port y   = {[round(y,3) for y in r['out_ports']]}")


if __name__ == "__main__":
    main()
