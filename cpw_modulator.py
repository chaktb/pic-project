"""LiNbO3 coplanar-waveguide (CPW) traveling-wave modulator calculator.

Port of the physics in public/cpw_modulator.html. Computes, over an RF sweep:
  * CPW RLGC line parameters and S11/S21 (electrical bandwidth),
  * RF attenuation and complex characteristic impedance,
  * microwave / optical velocity (effective-index) matching,
  * the traveling-wave electro-optic (EO) response, and
  * the half-wave voltage Vpi(f).

The optical waveguide index comes from Sellmeier dispersion of LiNbO3 plus a
Ti-indiffusion step (asymmetric slab, fundamental mode). Quasi-static CPW
parameters use the conformal-mapping (complete elliptic integral) formulas.

Usage:
    python cpw_modulator.py            # all defaults (matches the web form)
    python cpw_modulator.py --w 8 --s 16 --line_len 20 --pol TM
Run with -h for the full list of parameters (units shown in the help text).
"""

import argparse
import cmath
import math

# --------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------
MU0 = 4 * math.pi * 1e-7
EPS0 = 8.854187817e-12
C0 = 1 / math.sqrt(MU0 * EPS0)


# --------------------------------------------------------------------------
# Math utilities
# --------------------------------------------------------------------------
def linspace(a, b, n):
    if n <= 1:
        return [a]
    s = (b - a) / (n - 1)
    return [a + i * s for i in range(n)]


def interp(x, xs, ys):
    if len(xs) < 2:
        return ys[0] if ys else 0.0
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    lo, hi = 0, len(xs) - 1
    while hi - lo > 1:
        m = (lo + hi) >> 1
        if xs[m] <= x:
            lo = m
        else:
            hi = m
    t = (x - xs[lo]) / (xs[hi] - xs[lo])
    return ys[lo] + t * (ys[hi] - ys[lo])


def find_crossing(x, y, target):
    """First x where y crosses `target` (linear interpolation), else None."""
    for i in range(len(x) - 1):
        d1, d2 = y[i] - target, y[i + 1] - target
        if d1 == 0:
            return x[i]
        if d1 * d2 < 0:
            if y[i + 1] == y[i]:
                return x[i]
            return x[i] + (target - y[i]) * (x[i + 1] - x[i]) / (y[i + 1] - y[i])
    return None


def ellipk(m):
    """Complete elliptic integral K(m), m = k^2, via the AGM."""
    if m >= 1:
        return math.inf
    if m < 0:
        return math.nan
    a, b = 1.0, math.sqrt(1 - m)
    for _ in range(200):
        an, bn = 0.5 * (a + b), math.sqrt(a * b)
        if abs(an - bn) <= 1e-15 * an:
            a = an
            break
        a, b = an, bn
    return math.pi / (2 * a)


def brentq(f, xa, xb, xtol=1e-14, max_it=200):
    """Brent's root finder on [xa, xb]; returns None if not bracketed."""
    a, b = xa, xb
    fa, fb = f(a), f(b)
    if fa * fb > 0:
        return None
    if abs(fa) < abs(fb):
        a, b = b, a
        fa, fb = fb, fa
    c, fc, d, mflag = a, fa, a, True
    for _ in range(max_it):
        if abs(b - a) < xtol:
            return b
        if fb == 0:
            return b
        if fa != fc and fb != fc:
            s = (a * fb * fc / ((fa - fb) * (fa - fc))
                 + b * fa * fc / ((fb - fa) * (fb - fc))
                 + c * fa * fb / ((fc - fa) * (fc - fb)))
        else:
            s = b - fb * (b - a) / (fb - fa)
        bnd = 0.25 * (3 * a + b)
        c1 = (s - bnd) * (s - b) > 0
        c2 = mflag and abs(s - b) >= 0.5 * abs(b - c)
        c3 = (not mflag) and abs(s - b) >= 0.5 * abs(c - d)
        c4 = mflag and abs(b - c) < xtol
        c5 = (not mflag) and abs(c - d) < xtol
        if c1 or c2 or c3 or c4 or c5:
            s = 0.5 * (a + b)
            mflag = True
        else:
            mflag = False
        fs = f(s)
        d, c, fc = c, b, fb
        if fa * fs < 0:
            b, fb = s, fs
        else:
            a, fa = s, fs
        if abs(fa) < abs(fb):
            a, b = b, a
            fa, fb = fb, fa
    return b


# --------------------------------------------------------------------------
# Physics
# --------------------------------------------------------------------------
def sell_ln_no(l):
    L = l * l
    return math.sqrt(4.9048 + 0.11775 / (L - 0.0475) - 0.027169 * L)


def sell_ln_ne(l):
    L = l * l
    return math.sqrt(4.5820 + 0.09921 / (L - 0.04443) - 0.021950 * L)


def sell_sio2(l):
    L = l * l
    return math.sqrt(
        1
        + 0.6961663 * L / (L - 0.0684043 ** 2)
        + 0.4079426 * L / (L - 0.1162414 ** 2)
        + 0.8974794 * L / (L - 9.896161 ** 2)
    )


def slab_neff(nC, nS, nCl, d, lam):
    """Effective index of the fundamental mode of an asymmetric slab."""
    k0 = 2 * math.pi / lam
    n_lo = max(nS, nCl)
    if nC <= n_lo:
        return None
    dn2 = nC * nC - nS * nS
    if dn2 <= 0:
        return None
    V = k0 * d * math.sqrt(nC * nC - n_lo * n_lo)
    if V < 0.01:
        return None
    a = (nS * nS - nCl * nCl) / dn2
    Vs = k0 * d * math.sqrt(dn2)

    def eig(b):
        if b <= 0 or b >= 1:
            return 1e6
        sq = math.sqrt(1 - b)
        rhs = math.atan(math.sqrt(b / (1 - b)))
        a2 = (b + a) / (1 - b)
        if a2 < 0:
            return 1e6
        rhs += math.atan(math.sqrt(a2))
        return Vs * sq - rhs

    if eig(1e-12) * eig(1 - 1e-12) > 0:
        return nS
    b = brentq(eig, 1e-12, 1 - 1e-12, 1e-14)
    if b is None:
        return nS
    return math.sqrt(nS * nS + b * dn2)


def cpw_quasi_static(w, s, eff):
    """Quasi-static CPW: characteristic impedance and per-metre L, C."""
    k = w / (w + 2 * s)
    kp = math.sqrt(1 - k * k)
    Z0 = (30 * math.pi / math.sqrt(eff)) * (ellipk(kp * kp) / ellipk(k * k))
    vp = C0 / math.sqrt(eff)
    return Z0, 1 / (Z0 * vp), Z0 / vp  # Z0, C, L


def eff_eps(tox, s, eS, eL):
    """Effective RF permittivity (buffer-layer weighting)."""
    eta = 1 - math.exp(-2 * tox / s)
    return 0.5 * (1 + eta * eS + (1 - eta) * eL)


def r_pm(f, w, t, sig):
    """Per-metre conductor (skin-effect) resistance."""
    return 2.2 * math.sqrt(math.pi * f * MU0 / sig) / (w + 2 * t)


def g_pm(f, C, tS, tL, tox, s):
    """Per-metre dielectric (loss-tangent) conductance."""
    eta = 1 - math.exp(-2 * tox / s)
    return 2 * math.pi * f * C * (eta * tS + (1 - eta) * tL)


def tline_s(R, L, G, C, f, length, zref):
    """RLGC line -> (S11, S21, Zc, gamma) for a port impedance zref."""
    w = 2 * math.pi * f
    num = complex(R, w * L)
    den = complex(G, w * C)
    gamma = cmath.sqrt(num * den)
    zc = cmath.sqrt(num / den)
    gl = gamma * length
    A = cmath.cosh(gl)
    B = zc * cmath.sinh(gl)
    Cm = (1 / zc) * cmath.sinh(gl)
    D = cmath.cosh(gl)
    zr = complex(zref, 0)
    den2 = A + B / zr + Cm * zr + D
    s11 = (A + B / zr - Cm * zr - D) / den2
    s21 = 2 / den2
    return s11, s21, zc, gamma


def vpi_dc(lam, gap, n, r, G, Lel, pp):
    """DC half-wave voltage."""
    return lam * gap / ((2 if pp else 1) * n ** 3 * r * G * Lel)


def tw_eo(f, alpha, n_rf, n_opt, Lel):
    """Traveling-wave EO response magnitude (velocity / loss walk-off)."""
    d_beta = 2 * math.pi * f / C0 * (n_rf - n_opt)
    uL = complex(alpha, d_beta) * Lel
    if abs(uL) < 1e-12:
        return 1.0
    return abs((1 - cmath.exp(-uL)) / uL)


# --------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------
def simulate(p):
    f = linspace(p["f_start"], p["f_stop"], p["npts"])
    Lel = p["line_len"]

    eps_rf = eff_eps(p["tox"], p["s"], p["epsr_sio2"], p["epsr_ln"])
    Z0q, Cpm, Lpm = cpw_quasi_static(p["w"], p["s"], eps_rf)

    s11_db, s21_db = [], []
    alpha, zc_re, zc_im, n_rf, gammas = [], [], [], [], []
    for fi in f:
        R = r_pm(fi, p["w"], p["t_au"], p["sigma_au"])
        G = g_pm(fi, Cpm, p["tan_sio2"], p["tan_ln"], p["tox"], p["s"])
        s11, s21, zc, gamma = tline_s(R, Lpm, G, Cpm, fi, Lel, p["z0_port"])
        s11_db.append(20 * math.log10(max(abs(s11), 1e-15)))
        s21_db.append(20 * math.log10(max(abs(s21), 1e-15)))
        alpha.append(gamma.real * 8.686 / 100)            # Np/m -> dB/cm
        zc_re.append(zc.real)
        zc_im.append(zc.imag)
        n_rf.append(C0 * gamma.imag / (2 * math.pi * fi))
        gammas.append(gamma)

    f3 = find_crossing(f, s21_db, -3)
    f6 = find_crossing(f, s21_db, -6)
    n_rf_qs = math.sqrt(eps_rf)

    # Optical
    no_b = sell_ln_no(p["lam_um"])
    ne_b = sell_ln_ne(p["lam_um"])
    n_sio2 = sell_sio2(p["lam_um"])
    is_tm = p["pol"] == "TM"
    n_bulk = ne_b if is_tm else no_b
    pol_tag = "TM (ne)" if is_tm else "TE (no)"
    n_core = n_bulk + p["dn_ti"]
    n_opt = slab_neff(n_core, n_bulk, n_sio2, p["d_ti"], p["lam"])
    if n_opt is None:
        n_opt = n_bulk
    dl = 0.001
    sell = sell_ln_ne if is_tm else sell_ln_no
    dn_dl = (sell(p["lam_um"] + dl) - sell(p["lam_um"] - dl)) / (2 * dl)
    ng = n_bulk - p["lam_um"] * dn_dl

    # EO / Vpi
    r_eo = p["r33"] if is_tm else p["r13"]
    n_eo = ne_b if is_tm else no_b
    gap = p["s"]
    vpi_dc_v = vpi_dc(p["lam"], gap, n_eo, r_eo, p["gamma"], Lel, p["pp"])
    vpil_dc = vpi_dc_v * Lel * 100

    m_eo_db, vpif = [], []
    for i, fi in enumerate(f):
        m = tw_eo(fi, gammas[i].real, n_rf[i], n_opt, Lel)
        mc = max(m, 1e-15)
        m_eo_db.append(20 * math.log10(mc))
        vpif.append(vpi_dc_v / mc)
    f_eo3 = find_crossing(f, m_eo_db, -3)
    f_eo6 = find_crossing(f, m_eo_db, -6)

    return dict(
        f=f, eps_rf=eps_rf, Z0q=Z0q, Cpm=Cpm, Lpm=Lpm, n_rf_qs=n_rf_qs,
        s21_db=s21_db, alpha=alpha, zc_re=zc_re, n_rf=n_rf,
        f3=f3, f6=f6, no_b=no_b, ne_b=ne_b, n_sio2=n_sio2, n_core=n_core,
        n_opt=n_opt, ng=ng, pol_tag=pol_tag, is_tm=is_tm, r_eo=r_eo, n_eo=n_eo,
        gap=gap, Lel=Lel, vpi_dc=vpi_dc_v, vpil_dc=vpil_dc, vpif=vpif,
        f_eo3=f_eo3, f_eo6=f_eo6,
    )


def report(p, r):
    def fg(v):
        return "N/A" if v is None else f"{v / 1e9:.2f} GHz"

    f = r["f"]
    out = []
    out.append("== RF (CPW) ==========================")
    out.append(f"  eps_eff (RF)   = {r['eps_rf']:.4f}")
    out.append(f"  n_RF (q.s.)    = {r['n_rf_qs']:.4f}")
    out.append(f"  Z0 (q.s.)      = {r['Z0q']:.2f} ohm")
    out.append(f"  L              = {r['Lpm'] * 1e6:.4f} uH/m")
    out.append(f"  C              = {r['Cpm'] * 1e12:.4f} pF/m")
    out.append(f"  v_ph (RF)      = {C0 / r['n_rf_qs'] / 1e8:.4f} x10^8 m/s")
    out.append("")
    out.append(f"  -3 dB freq     = {fg(r['f3'])}")
    out.append(f"  -6 dB freq     = {fg(r['f6'])}")
    out.append(f"  alpha @ 10 GHz = {interp(10e9, f, r['alpha']):.4f} dB/cm")
    out.append(f"  alpha @ 40 GHz = {interp(40e9, f, r['alpha']):.4f} dB/cm")
    out.append(f"  Re(Zc)@10G     = {interp(10e9, f, r['zc_re']):.2f} ohm")
    out.append("")
    out.append("== Optical (LiNbO3) ==================")
    out.append(f"  lambda         = {p['lam_um']} um")
    out.append(f"  Pol.           = {r['pol_tag']}")
    out.append(f"  no (bulk)      = {r['no_b']:.5f}")
    out.append(f"  ne (bulk)      = {r['ne_b']:.5f}")
    out.append(f"  n_SiO2         = {r['n_sio2']:.5f}")
    out.append(f"  n_core (+dn)   = {r['n_core']:.5f}")
    out.append(f"  n_eff (WG)     = {r['n_opt']:.5f}")
    out.append(f"  n_g (group)    = {r['ng']:.5f}")
    out.append("")
    out.append("== Velocity Matching =================")
    dn_mid = interp(0.5 * (p["f_start"] + p["f_stop"]), f, r["n_rf"]) - r["n_opt"]
    out.append(f"  dn (mid-freq)  = {dn_mid:+.5f}")
    dn_qs = r["n_rf_qs"] - r["n_opt"]
    out.append(f"  n_RF - n_opt   = {dn_qs:+.5f}  (q.s.)")
    out.append("")
    out.append("== Vpi / EO ==========================")
    out.append(f"  EO coeff       = {'r33' if r['is_tm'] else 'r13'} = {r['r_eo'] * 1e12:.1f} pm/V")
    out.append(f"  n (EO)         = {r['n_eo']:.5f}")
    out.append(f"  Gamma (overlap)= {p['gamma']:.3f}")
    out.append(f"  MZM type       = {'Push-pull' if p['pp'] else 'Single-arm'}")
    out.append(f"  Electrode gap  = {r['gap'] * 1e6:.1f} um")
    out.append(f"  Electrode L    = {r['Lel'] * 1e3:.1f} mm")
    out.append("")
    out.append(f"  Vpi (DC)       = {r['vpi_dc']:.3f} V")
    out.append(f"  Vpi.L (DC)     = {r['vpil_dc']:.3f} V.cm")
    out.append("")
    out.append(f"  Vpi @ 10 GHz   = {interp(10e9, f, r['vpif']):.3f} V")
    out.append(f"  Vpi @ 40 GHz   = {interp(40e9, f, r['vpif']):.3f} V")
    out.append("")
    out.append(f"  EO BW (-3dB el)= {fg(r['f_eo3'])}")
    out.append(f"  EO BW (-6dB op)= {fg(r['f_eo6'])}")
    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser(
        description="LiNbO3 CPW traveling-wave modulator calculator")
    # CPW geometry
    ap.add_argument("--w", type=float, default=8.0, help="signal width [um]")
    ap.add_argument("--s", type=float, default=16.0, help="gap [um]")
    ap.add_argument("--t_au", type=float, default=20.0, help="metal thickness [um]")
    ap.add_argument("--tox", type=float, default=1.5, help="SiO2 thickness [um]")
    ap.add_argument("--line_len", type=float, default=20.0, help="electrode length [mm]")
    # RF materials
    ap.add_argument("--sigma_au", type=float, default=4.1e7, help="sigma_Au [S/m]")
    ap.add_argument("--epsr_sio2", type=float, default=3.9)
    ap.add_argument("--epsr_ln", type=float, default=28.0)
    ap.add_argument("--tan_sio2", type=float, default=0.0005)
    ap.add_argument("--tan_ln", type=float, default=0.003)
    # Optical
    ap.add_argument("--lam", type=float, default=1.55, help="wavelength [um]")
    ap.add_argument("--dn_ti", type=float, default=0.008, help="Ti index step")
    ap.add_argument("--d_ti", type=float, default=4.0, help="Ti depth [um]")
    ap.add_argument("--pol", choices=["TM", "TE"], default="TM")
    # EO / Vpi
    ap.add_argument("--r33", type=float, default=30.8, help="r33 [pm/V]")
    ap.add_argument("--r13", type=float, default=8.6, help="r13 [pm/V]")
    ap.add_argument("--gamma", type=float, default=0.5, help="overlap factor")
    ap.add_argument("--single-arm", dest="pp", action="store_false",
                    help="single-arm (default is push-pull)")
    ap.set_defaults(pp=True)
    # Frequency
    ap.add_argument("--f_start", type=float, default=1.0, help="start [GHz]")
    ap.add_argument("--f_stop", type=float, default=67.0, help="stop [GHz]")
    ap.add_argument("--npts", type=int, default=800)
    ap.add_argument("--z0_port", type=float, default=50.0, help="port Z0 [ohm]")
    a = ap.parse_args()

    p = dict(
        w=a.w * 1e-6, s=a.s * 1e-6, t_au=a.t_au * 1e-6, tox=a.tox * 1e-6,
        line_len=a.line_len * 1e-3,
        sigma_au=a.sigma_au, epsr_sio2=a.epsr_sio2, epsr_ln=a.epsr_ln,
        tan_sio2=a.tan_sio2, tan_ln=a.tan_ln,
        lam_um=a.lam, lam=a.lam * 1e-6, dn_ti=a.dn_ti, d_ti=a.d_ti * 1e-6,
        pol=a.pol, r33=a.r33 * 1e-12, r13=a.r13 * 1e-12, gamma=a.gamma, pp=a.pp,
        f_start=a.f_start * 1e9, f_stop=a.f_stop * 1e9,
        npts=max(2, a.npts), z0_port=a.z0_port,
    )
    print(report(p, simulate(p)))


if __name__ == "__main__":
    main()
