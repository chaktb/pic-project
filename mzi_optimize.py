"""MZI optimizer.

Given a list of target points (lambda_i, P_bar_target_i, weight_i), search for
DC parameters (gap1, L1, gap2, L2, Delta_L) that minimize the weighted sum of
squared errors. Uses multi-start Nelder-Mead.

Usage:
    python mzi_optimize.py <W> <H> <delta_pct> <targets_csv> [--starts N]

targets_csv: comma-separated 'lambda:P_bar[:weight]' (e.g., '0.98:1.0,1.55:0.0')

Example (980/1550 WDM):
    python mzi_optimize.py 5.0 5.0 0.75 0.98:1.0,1.55:0.0
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
    return 2 * h2 * gx * math.exp(-gx * gap) / (beta * (h2 + gx*gx) * (1 + gx * a)), n_eff


def p_bar(wl, W, H, delta, gap1, L1, gap2, L2, dL):
    kn1 = kappa_neff(wl, W, H, delta, gap1)
    kn2 = kappa_neff(wl, W, H, delta, gap2)
    if kn1 is None or kn2 is None:
        return None
    k1, n_eff = kn1
    k2, _     = kn2
    kL1 = k1 * L1
    kL2 = k2 * L2
    dphi = (2 * math.pi * n_eff / wl) * dL
    a_ = math.cos(kL1) * math.cos(kL2)
    b_ = math.sin(kL1) * math.sin(kL2)
    return a_*a_ + b_*b_ - 2*a_*b_*math.cos(dphi)


def make_cost(W, H, delta, targets, bounds):
    """targets = [(wl, p_target, weight), ...].  bounds = [(lo, hi), ...] for 5 params."""
    def clip(x):
        return [max(lo, min(hi, xi)) for xi, (lo, hi) in zip(x, bounds)]

    def cost(x):
        g1, L1, g2, L2, dL = clip(x)
        s = 0.0
        for wl, t, w in targets:
            p = p_bar(wl, W, H, delta, g1, L1, g2, L2, dL)
            if p is None:
                return 1e6
            s += w * (p - t) ** 2
        # Soft penalty for hitting bounds
        for xi, xc, (lo, hi) in zip(x, clip(x), bounds):
            s += 0.001 * ((xi - xc) / max(hi - lo, 1e-12)) ** 2
        return s
    return cost


def nelder_mead(f, x0, max_iter=200, tol=1e-8):
    n = len(x0)
    simplex = [list(x0)]
    for i in range(n):
        xi = list(x0)
        step = abs(xi[i]) * 0.05 if abs(xi[i]) > 1e-6 else 0.001
        xi[i] = xi[i] + step
        simplex.append(xi)
    fvals = [f(p) for p in simplex]
    for _ in range(max_iter):
        order = sorted(range(n+1), key=lambda i: fvals[i])
        simplex = [simplex[i] for i in order]
        fvals = [fvals[i] for i in order]
        if fvals[n] - fvals[0] < tol:
            break
        c = [sum(simplex[i][j] for i in range(n)) / n for j in range(n)]
        # Reflect
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
        xc = [c[j] + 0.5 * (simplex[n][j] - c[j]) for j in range(n)]
        fc = f(xc)
        if fc < fvals[n]:
            simplex[n], fvals[n] = xc, fc
            continue
        # Shrink
        for i in range(1, n+1):
            simplex[i] = [simplex[0][j] + 0.5 * (simplex[i][j] - simplex[0][j]) for j in range(n)]
            fvals[i] = f(simplex[i])
    order = sorted(range(n+1), key=lambda i: fvals[i])
    return simplex[order[0]], fvals[order[0]]


def optimize(W, H, delta, targets, bounds=None, n_starts=30, seed=1):
    if bounds is None:
        bounds = [(1.0, 15.0), (10.0, 5000.0), (1.0, 15.0), (10.0, 5000.0), (0.1, 200.0)]
    rnd = random.Random(seed)
    cost = make_cost(W, H, delta, targets, bounds)
    best_x, best_c = None, float("inf")
    for _ in range(n_starts):
        x0 = [lo + rnd.random() * (hi - lo) for (lo, hi) in bounds]
        x, c = nelder_mead(cost, x0, max_iter=300)
        x = [max(lo, min(hi, xi)) for xi, (lo, hi) in zip(x, bounds)]
        c2 = cost(x)
        if c2 < best_c:
            best_c, best_x = c2, x
    return best_x, best_c


def parse_targets(s):
    out = []
    for part in s.split(","):
        part = part.strip()
        if not part:
            continue
        bits = part.split(":")
        wl = float(bits[0])
        t = float(bits[1])
        w = float(bits[2]) if len(bits) >= 3 else 1.0
        out.append((wl, t, w))
    return out


def main():
    if len(sys.argv) < 5:
        print("Usage: python mzi_optimize.py <W> <H> <delta_pct> <targets> [--starts N]")
        print("Example: python mzi_optimize.py 5.0 5.0 0.75 0.98:1.0,1.55:0.0")
        sys.exit(1)
    W = float(sys.argv[1])
    H = float(sys.argv[2])
    dl = float(sys.argv[3])
    targets = parse_targets(sys.argv[4])
    n_starts = 30
    for i, a in enumerate(sys.argv):
        if a == "--starts" and i+1 < len(sys.argv):
            n_starts = int(sys.argv[i+1])

    print(f"# optimizing for {len(targets)} targets, n_starts={n_starts}")
    x, c = optimize(W, H, dl, targets, n_starts=n_starts)
    g1, L1, g2, L2, dL = x
    print(f"best cost = {c:.6e}")
    print(f"DC1: gap = {g1:.4f} um   L = {L1:.4f} um")
    print(f"DC2: gap = {g2:.4f} um   L = {L2:.4f} um")
    print(f"Delta_L  = {dL:.4f} um")
    print()
    print("verification (target -> achieved):")
    for wl, t, w in targets:
        p = p_bar(wl, W, H, dl, g1, L1, g2, L2, dL)
        print(f"  lambda={wl:.4f} um  P_bar target={t:.3f}  achieved={p:.6f}  P_cross={1-p:.6f}")


if __name__ == "__main__":
    main()
