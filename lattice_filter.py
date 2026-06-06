"""Lattice-form bandpass filter optimizer.

An N-stage lattice filter is (N+1) directional couplers separated by N
equal-length delay lines (each with path difference dL on the bottom arm).
A 1-stage lattice is just an MZI; higher N gives flatter passband and
sharper skirts.

Usage:
    python lattice_filter.py <W> <H> <delta_pct> <wl_c_um> <passband_um> <N_stages>
Example (1550 nm center, 0.2 nm passband, 3 stages):
    python lattice_filter.py 5.0 5.0 0.75 1.55 0.0002 3
"""

import math
import random
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
    k = 2 * h2 * gx * math.exp(-gx * gap) / (beta * (h2 + gx*gx) * (1 + gx * a))
    return k, n_eff


def lattice_pbar(wl, W, H, delta, gap, Ls, dL):
    """Bar power for an N-stage lattice (len(Ls) = N+1)."""
    kn = kappa_neff(wl, W, H, delta, gap)
    if kn is None:
        return None
    kappa, n_eff = kn
    dphi = (2 * math.pi * n_eff / wl) * dL
    cp, sp = math.cos(dphi), math.sin(dphi)

    uR, uI, vR, vI = 1.0, 0.0, 0.0, 0.0
    N = len(Ls) - 1
    for k in range(N + 1):
        theta = kappa * Ls[k]
        c, s = math.cos(theta), math.sin(theta)
        nuR = c*uR + s*vI
        nuI = c*uI - s*vR
        nvR = c*vR + s*uI
        nvI = c*vI - s*uR
        uR, uI, vR, vI = nuR, nuI, nvR, nvI
        if k < N:
            nvR2 = vR*cp + vI*sp
            nvI2 = -vR*sp + vI*cp
            vR, vI = nvR2, nvI2
    return uR*uR + uI*uI


def bandpass_targets(wl_c, passband):
    """Generate target points for a bandpass filter centered at wl_c (μm)
    with FWHM passband Δλ (μm)."""
    return [
        (wl_c,                 1.0, 5.0),
        (wl_c - passband/2,    0.5, 1.0),
        (wl_c + passband/2,    0.5, 1.0),
        (wl_c - passband,      0.05, 2.0),
        (wl_c + passband,      0.05, 2.0),
        (wl_c - 2*passband,    0.0, 1.0),
        (wl_c + 2*passband,    0.0, 1.0),
    ]


def make_cost(W, H, delta, gap, targets, n_stages, bounds):
    def clip(x):
        return [max(lo, min(hi, xi)) for xi, (lo, hi) in zip(x, bounds)]

    def cost(x):
        xc = clip(x)
        Ls = xc[:n_stages + 1]
        dL = xc[n_stages + 1]
        s = 0.0
        for wl, t, w in targets:
            p = lattice_pbar(wl, W, H, delta, gap, Ls, dL)
            if p is None:
                return 1e6
            s += w * (p - t) ** 2
        for xi, xci, (lo, hi) in zip(x, xc, bounds):
            s += 0.001 * ((xi - xci) / max(hi - lo, 1e-12)) ** 2
        return s
    return cost


def nelder_mead(f, x0, max_iter=300, tol=1e-9):
    n = len(x0)
    simplex = [list(x0)]
    for i in range(n):
        xi = list(x0)
        step = abs(xi[i]) * 0.05 if abs(xi[i]) > 1e-6 else 0.001
        xi[i] += step
        simplex.append(xi)
    fvals = [f(p) for p in simplex]
    for _ in range(max_iter):
        order = sorted(range(n+1), key=lambda i: fvals[i])
        simplex = [simplex[i] for i in order]
        fvals = [fvals[i] for i in order]
        if fvals[n] - fvals[0] < tol:
            break
        c = [sum(simplex[i][j] for i in range(n)) / n for j in range(n)]
        xr = [c[j] + (c[j] - simplex[n][j]) for j in range(n)]
        fr = f(xr)
        if fvals[0] <= fr < fvals[n-1]:
            simplex[n], fvals[n] = xr, fr
            continue
        if fr < fvals[0]:
            xe = [c[j] + 2 * (xr[j] - c[j]) for j in range(n)]
            fe = f(xe)
            if fe < fr:
                simplex[n], fvals[n] = xe, fe
            else:
                simplex[n], fvals[n] = xr, fr
            continue
        xc2 = [c[j] + 0.5 * (simplex[n][j] - c[j]) for j in range(n)]
        fc = f(xc2)
        if fc < fvals[n]:
            simplex[n], fvals[n] = xc2, fc
            continue
        for i in range(1, n+1):
            simplex[i] = [simplex[0][j] + 0.5 * (simplex[i][j] - simplex[0][j]) for j in range(n)]
            fvals[i] = f(simplex[i])
    order = sorted(range(n+1), key=lambda i: fvals[i])
    return simplex[order[0]], fvals[order[0]]


def optimize(W, H, delta, gap, targets, n_stages, bounds=None,
             n_starts=80, seed=1):
    if bounds is None:
        bounds = [(10.0, 5000.0)] * (n_stages + 1) + [(50.0, 5000.0)]
    rnd = random.Random(seed)
    cost = make_cost(W, H, delta, gap, targets, n_stages, bounds)
    best_x, best_c = None, float("inf")
    for _ in range(n_starts):
        x0 = [lo + rnd.random() * (hi - lo) for (lo, hi) in bounds]
        x, _ = nelder_mead(cost, x0, max_iter=400)
        x = [max(lo, min(hi, xi)) for xi, (lo, hi) in zip(x, bounds)]
        c = cost(x)
        if c < best_c:
            best_c, best_x = c, x
    return best_x, best_c


def main():
    if len(sys.argv) < 7:
        print("Usage: python lattice_filter.py <W> <H> <delta_pct> <wl_c_um> <passband_um> <N_stages>")
        print("Example: python lattice_filter.py 5.0 5.0 0.75 1.55 0.0002 3")
        sys.exit(1)
    W   = float(sys.argv[1])
    H   = float(sys.argv[2])
    dl  = float(sys.argv[3])
    wlc = float(sys.argv[4])
    pb  = float(sys.argv[5])
    N   = int(sys.argv[6])
    gap = 5.0

    targets = bandpass_targets(wlc, pb)
    x, c = optimize(W, H, dl, gap, targets, N, n_starts=80, seed=1)
    Ls = x[:N+1]
    dL = x[N+1]

    print(f"# N stages = {N},  center = {wlc} um,  passband = {pb*1000:.3f} nm")
    print(f"best cost = {c:.6e}")
    for i, L in enumerate(Ls):
        print(f"  DC{i}: gap = {gap:.3f} um   L = {L:.4f} um")
    print(f"  dL  = {dL:.4f} um  (per stage)")
    print()
    print("spectrum near center:")
    for off in (-2, -1, -0.5, 0, 0.5, 1, 2):
        wl = wlc + off * pb
        p = lattice_pbar(wl, W, H, dl, gap, Ls, dL)
        print(f"  lambda={wl*1000:.3f} nm  ({off:+.1f}·Δλ)  P_bar={p:.6f}")


if __name__ == "__main__":
    main()
