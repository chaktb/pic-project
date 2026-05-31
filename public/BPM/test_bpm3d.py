"""Headless numerical verification for BPM3DSolver (no GUI, no GDS)."""

import numpy as np
from bpm3d import BPM3DSolver, VerticalStack, StackLayer, default_silica_stack


def make_solver(n_ref, dx=0.1, dy=0.1, dz=0.2, lam=1.55, pml=2.0):
    return BPM3DSolver(wavelength_um=lam, n_ref=n_ref, dx_um=dx, dy_um=dy,
                       dz_um=dz, pml_um=pml, pml_sigma=4.0, pml_m=3.0)


def beam_width(profile, coord):
    """1/e^2 intensity half-width from second moment: w = 2*sqrt(<c^2>)."""
    p = profile / np.trapz(profile, coord)
    c0 = np.trapz(coord * p, coord)
    var = np.trapz((coord - c0) ** 2 * p, coord)
    return 2.0 * np.sqrt(max(var, 0.0)), c0


# ---------------------------------------------------------------------------
def test_free_space_diffraction():
    """Uniform medium: Gaussian beam must diffract per analytic Gaussian optics."""
    n = 1.444
    s = make_solver(n_ref=n, dx=0.1, dy=0.1, dz=0.1, pml=3.0)
    x = s.make_symmetric_grid(-20, 20, s.dx_um)
    y = s.make_symmetric_grid(-20, 20, s.dy_um)
    w0 = 2.0  # amplitude 1/e half-width
    field0 = s.gaussian_2d(x, y, 0, 0, w0, w0)

    Lz = 30.0
    z = np.arange(0, Lz + 0.5 * s.dz_um, s.dz_um)
    idx = np.full((len(x), len(y)), n)
    res = s.run(x, y, z, lambda iz: idx, field0, save_every=8, progress_cb=None)

    # For intensity exp(-2x^2/w0^2) the second-moment width 2*sqrt(var) == w0,
    # so the analytic width is simply w0 * sqrt(1+(z/zR)^2).
    zR = np.pi * w0 ** 2 * n / s.wavelength_um
    Wend_analytic = w0 * np.sqrt(1 + (Lz / zR) ** 2)

    cut = res.psi_out
    px = np.trapz(np.abs(cut) ** 2, y, axis=1)
    Wend_sim, _ = beam_width(px, x)

    err = abs(Wend_sim - Wend_analytic) / Wend_analytic
    p0, p1 = res.total_power[0], res.total_power[-1]
    print(f"[free-space] zR={zR:.2f}um  W_analytic={Wend_analytic:.3f}  "
          f"W_sim={Wend_sim:.3f}  err={err*100:.1f}%  "
          f"power {p0:.4f}->{p1:.4f}")
    assert err < 0.06, f"diffraction width off by {err*100:.1f}%"
    assert abs(p1 - p0) / p0 < 0.05, "free-space power not conserved"
    print("  PASS")


def test_straight_waveguide():
    """Buried square core: a launched beam stays guided, power ~ conserved."""
    n_clad, n_core = 1.444, 1.4549   # ~0.75% contrast
    wx_core, wy_core = 6.0, 6.0
    # crude guided n_ref between clad and core
    n_ref = 0.5 * (n_core + n_clad) + 0.4 * (n_core - n_clad)
    s = make_solver(n_ref=n_ref, dx=0.15, dy=0.15, dz=0.5, pml=3.0)
    x = s.make_symmetric_grid(-18, 18, s.dx_um)
    y = s.make_symmetric_grid(-18, 18, s.dy_um)

    core = (np.abs(x)[:, None] <= wx_core / 2) & (np.abs(y)[None, :] <= wy_core / 2)
    idx = np.where(core, n_core, n_clad)

    field0 = s.gaussian_2d(x, y, 0, 0, wx_core / 2.2, wy_core / 2.2)
    Lz = 200.0
    z = np.arange(0, Lz + 0.5 * s.dz_um, s.dz_um)
    res = s.run(x, y, z, lambda iz: idx, field0, save_every=10)

    # confined power in a window around the core, start vs end
    win_x = np.abs(x) <= wx_core
    win_y = np.abs(y) <= wy_core
    def conf(fld):
        sub = fld[np.ix_(win_x, win_y)]
        return float(np.trapz(np.trapz(np.abs(sub) ** 2, y[win_y], axis=1), x[win_x]))
    c0 = conf(field0)
    c1 = conf(res.psi_out)
    p0, p1 = res.total_power[0], res.total_power[-1]
    wsim, _ = beam_width(np.trapz(np.abs(res.psi_out) ** 2, y, axis=1), x)
    print(f"[straight wg] total power {p0:.4f}->{p1:.4f}  "
          f"confined {c0:.4f}->{c1:.4f} ({c1/c0*100:.1f}%)  out width={wsim:.2f}um")
    assert abs(p1 - p0) / p0 < 0.08, "total power drifted too much"
    assert c1 / c0 > 0.80, "lost guiding — beam not confined"
    assert wsim < wx_core * 1.6, "beam spread far beyond core (not guided)"
    print("  PASS")


def test_pml_absorption():
    """A beam aimed at the wall must be absorbed, not reflected back."""
    n = 1.444
    s = make_solver(n_ref=n, dx=0.1, dy=0.1, dz=0.1, pml=3.0)
    x = s.make_symmetric_grid(-15, 15, s.dx_um)
    y = s.make_symmetric_grid(-15, 15, s.dy_um)
    # tilt the beam in +x by adding linear phase
    w0 = 1.5
    field0 = s.gaussian_2d(x, y, 0, 0, w0, w0)
    theta = np.deg2rad(10.0)
    field0 = field0 * np.exp(1j * s.beta_ref * np.sin(theta) * x[:, None])
    field0 = s.normalize(field0, x, y)

    Lz = 100.0  # at 10deg the beam centre walks ~18um in +x, well past the +x wall
    z = np.arange(0, Lz + 0.5 * s.dz_um, s.dz_um)
    idx = np.full((len(x), len(y)), n)
    res = s.run(x, y, z, lambda iz: idx, field0, save_every=8)

    p0, p1 = res.total_power[0], res.total_power[-1]
    # most power should have left through the PML wall
    print(f"[pml] power {p0:.4f}->{p1:.4f}  (absorbed {(1-p1/p0)*100:.0f}%)")
    assert p1 / p0 < 0.5, "beam not absorbed at wall — PML too weak/reflecting"
    # no large reflected lobe near the opposite (-x) wall
    px = np.trapz(np.abs(res.psi_out) ** 2, y, axis=1)
    left_lobe = px[x < x[0] + 4.0].sum()
    assert left_lobe < 0.02 * px.sum() + 1e-12, "reflection lobe at far wall"
    print("  PASS")


def analytic_slab_neff(a, n1, n2, k0, pol):
    """Fundamental even mode n_eff of a symmetric slab (half-thickness a).
    TE: u tan u = sqrt(V^2-u^2);  TM: u tan u = (n1/n2)^2 sqrt(V^2-u^2)."""
    V = k0 * a * np.sqrt(n1 ** 2 - n2 ** 2)
    fac = (n1 ** 2 / n2 ** 2) if pol == "TM" else 1.0

    def f(u):
        return u * np.tan(u) - fac * np.sqrt(max(V ** 2 - u ** 2, 0.0))

    lo, hi = 1e-7, min(V, np.pi / 2) - 1e-7
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(lo) * f(mid) <= 0:
            hi = mid
        else:
            lo = mid
    u = 0.5 * (lo + hi)
    beta = np.sqrt((n1 * k0) ** 2 - (u / a) ** 2)
    return beta / k0


def test_mode_slab_te_tm():
    """Mode solver on a y-slab must reproduce analytic slab TE0/TM0 n_eff and
    the polarization splitting (scalar == TE; TM lower than TE)."""
    n_clad, n_core, t = 1.444, 1.60, 4.0     # high contrast -> visible TE/TM split
    lam = 1.55
    dx, dy = 0.25, 0.04                      # fine in y (sharp TM interface), coarse in x
    k0 = 2 * np.pi / lam
    te0 = analytic_slab_neff(t / 2, n_core, n_clad, k0, "TE")
    tm0 = analytic_slab_neff(t / 2, n_core, n_clad, k0, "TM")

    # wide-x window (≈ free in x), slab varies only in y; NO PML so the x box
    # mode is a clean low-k cosine shared by every polarization.
    x = BPM3DSolver.make_symmetric_grid(-16, 16, dx)
    y = BPM3DSolver.make_symmetric_grid(-7, 7, dy)
    band = np.abs(y) <= t / 2
    n2 = np.where(band[None, :], n_core ** 2, n_clad ** 2) * np.ones((len(x), 1))

    def solve(pol):
        s = BPM3DSolver(wavelength_um=lam, n_ref=0.5 * (n_core + n_clad),
                        dx_um=dx, dy_um=dy, dz_um=0.2, pml_um=0.0, pol=pol)
        _, neff = s.solve_mode(x, y, n2, 0.0, 0.0, 10.0, t / 2.2, max_iter=600, tol=1e-8)
        return neff

    ns, nte, ntm = solve("scalar"), solve("TE"), solve("TM")
    print(f"[mode slab] analytic TE0={te0:.5f} TM0={tm0:.5f} | "
          f"solver scalar={ns:.5f} TE={nte:.5f} TM={ntm:.5f}")
    # the shared x box-mode adds an identical tiny offset to all three, so the
    # polarization *splitting* is the sharp check; absolute within a few e-3.
    assert abs(nte - te0) < 3e-3, f"TE n_eff off: {nte} vs {te0}"
    assert abs(ns - te0) < 3e-3, f"scalar n_eff should match TE0: {ns} vs {te0}"
    assert abs(ntm - tm0) < 3e-3, f"TM n_eff off: {ntm} vs {tm0}"
    assert abs((nte - ntm) - (te0 - tm0)) < 1e-3, "TE-TM splitting wrong"
    assert ntm < nte, "TM must be less confined than TE here"
    print("  PASS")


def test_eigenmode_launch():
    """A solved eigenmode launched into a straight waveguide should keep its
    shape (near-unit self-overlap) — far less breathing than a Gaussian."""
    n_clad, n_core, w = 1.444, 1.4549, 6.0
    lam = 1.55
    dx = dy = 0.15
    n_ref = 0.5 * (n_core + n_clad) + 0.4 * (n_core - n_clad)
    s = BPM3DSolver(wavelength_um=lam, n_ref=n_ref, dx_um=dx, dy_um=dy, dz_um=0.5,
                    pml_um=3.0, pol="scalar")
    x = s.make_symmetric_grid(-15, 15, dx)
    y = s.make_symmetric_grid(-15, 15, dy)
    core = (np.abs(x)[:, None] <= w / 2) & (np.abs(y)[None, :] <= w / 2)
    n2 = np.where(core, n_core ** 2, n_clad ** 2)
    idx = np.sqrt(n2)

    mode, neff = s.solve_mode(x, y, n2, 0, 0, w / 2.2, w / 2.2, max_iter=600)
    assert n_clad < neff < n_core, f"mode n_eff {neff} out of (n_clad, n_core)"

    Lz = 120.0
    z = np.arange(0, Lz + 0.5 * s.dz_um, s.dz_um)
    res_mode = s.run(x, y, z, lambda iz: idx, mode, save_every=10)
    gauss = s.gaussian_2d(x, y, 0, 0, w / 2.2, w / 2.2)
    res_g = s.run(x, y, z, lambda iz: idx, gauss, save_every=10)

    ov_mode = s.overlap_power(res_mode.psi_out, mode, x, y)
    ov_g = s.overlap_power(res_g.psi_out, gauss, x, y)
    print(f"[eigenmode] n_eff={neff:.5f}  self-overlap mode={ov_mode*100:.1f}%  "
          f"gaussian={ov_g*100:.1f}%  power={res_mode.total_power[-1]:.3f}")
    assert ov_mode > 0.97, f"eigenmode should stay coherent, got {ov_mode}"
    assert ov_mode > ov_g, "eigenmode must breathe less than the Gaussian"
    print("  PASS")


def _comp_power(fld, x, y):
    return float(np.trapz(np.trapz(np.abs(fld) ** 2, y, axis=1), x))


def test_fullvec_separable_no_coupling():
    """A separable (x-uniform) y-slab has ∂ln n²/∂x = 0 everywhere, so a pure-Ex
    launch must generate exactly zero Ey under full-vector propagation."""
    n_clad, n_core, t = 1.444, 1.60, 4.0
    dx = dy = 0.1
    s = BPM3DSolver(wavelength_um=1.55, n_ref=0.5 * (n_core + n_clad),
                    dx_um=dx, dy_um=dy, dz_um=0.4, pml_um=2.0, pol="full")
    x = s.make_symmetric_grid(-10, 10, dx)
    y = s.make_symmetric_grid(-8, 8, dy)
    band = np.abs(y) <= t / 2
    n2 = np.where(band[None, :], n_core ** 2, n_clad ** 2) * np.ones((len(x), 1))
    idx = np.sqrt(n2)

    ex0 = s.gaussian_2d(x, y, 0, 0, 2.0, t / 2.2)
    ey0 = np.zeros_like(ex0)
    z = np.arange(0, 40 + 0.2, s.dz_um)
    res = s.run_fullvec(x, y, z, lambda iz: idx, ex0, ey0, save_every=10)
    p_ey = _comp_power(res.psi_out_minor, x, y)
    p_tot = p_ey + _comp_power(res.psi_out, x, y)
    print(f"[fullvec separable] Ey/total = {p_ey / p_tot:.2e} (should be ~0)")
    assert p_ey / p_tot < 1e-10, "separable structure must not couple Ex->Ey"
    print("  PASS")


def test_fullvec_power_conservation():
    """Total power |Ex|²+|Ey|² conserved for a guided mixed-polarization launch."""
    n_clad, n_core, w = 1.444, 1.4549, 6.0
    dx = dy = 0.15
    n_ref = 0.5 * (n_core + n_clad) + 0.4 * (n_core - n_clad)
    s = BPM3DSolver(wavelength_um=1.55, n_ref=n_ref, dx_um=dx, dy_um=dy, dz_um=0.5,
                    pml_um=3.0, pol="full")
    x = s.make_symmetric_grid(-15, 15, dx)
    y = s.make_symmetric_grid(-15, 15, dy)
    core = (np.abs(x)[:, None] <= w / 2) & (np.abs(y)[None, :] <= w / 2)
    idx = np.where(core, n_core, n_clad)
    ex0 = s.gaussian_2d(x, y, 0, 0, w / 2.2, w / 2.2)
    ey0 = s.gaussian_2d(x, y, 0, 0, w / 2.2, w / 2.2)
    z = np.arange(0, 120 + 0.25, s.dz_um)
    res = s.run_fullvec(x, y, z, lambda iz: idx, ex0, ey0, save_every=10)
    p0, p1 = res.total_power[0], res.total_power[-1]
    print(f"[fullvec power] total {p0:.4f}->{p1:.4f}  conversion={res.conversion*100:.2f}%")
    assert abs(p1 - p0) / p0 < 0.08, "full-vector total power not conserved"
    print("  PASS")


def test_fullvec_rotated_conversion():
    """A 45°-rotated (diamond) core converts Ex->Ey; an axis-aligned square
    barely does. Confirms the cross-coupling term is active and physical."""
    n_clad, n_core, a = 1.444, 1.70, 1.6       # high contrast, small core
    dx = dy = 0.08
    n_ref = 0.5 * (n_core + n_clad) + 0.4 * (n_core - n_clad)

    def run_core(core_mask, Lz=45.0):
        s = BPM3DSolver(wavelength_um=1.55, n_ref=n_ref, dx_um=dx, dy_um=dy, dz_um=0.3,
                        pml_um=2.0, pol="full")
        idx = np.where(core_mask, n_core, n_clad)
        ex0 = s.gaussian_2d(x, y, 0, 0, a, a)   # pure Ex
        ey0 = np.zeros_like(ex0)
        z = np.arange(0, Lz + 0.15, s.dz_um)
        res = s.run_fullvec(x, y, z, lambda iz: idx, ex0, ey0, save_every=20)
        return res.conversion

    sg = BPM3DSolver(wavelength_um=1.55, n_ref=n_ref, dx_um=dx, dy_um=dy, dz_um=0.3, pml_um=2.0)
    x = sg.make_symmetric_grid(-7, 7, dx)
    y = sg.make_symmetric_grid(-7, 7, dy)
    aligned = (np.abs(x)[:, None] <= a) & (np.abs(y)[None, :] <= a)
    xr = (x[:, None] + y[None, :]) / np.sqrt(2)
    yr = (-x[:, None] + y[None, :]) / np.sqrt(2)
    diamond = (np.abs(xr) <= a) & (np.abs(yr) <= a)

    c_align = run_core(aligned)
    c_diam = run_core(diamond)
    print(f"[fullvec conversion] aligned Ex->Ey={c_align*100:.3f}%  "
          f"diamond(45°)={c_diam*100:.3f}%  ratio={c_diam/max(c_align,1e-9):.1f}x")
    assert c_diam > 2e-3, f"45° core should convert measurably, got {c_diam}"
    assert c_diam > 4 * c_align, "rotated core must convert much more than aligned"
    print("  PASS")


def test_fullvec_hybrid_mode():
    """The full-vector hybrid eigenmode (solve_mode_fullvec) must be stationary:
    launched into a straight guide it keeps its (Ex,Ey) shape (near-unit vector
    self-overlap), beating a non-eigen Gaussian launch."""
    n_clad, n_core, wx, wy = 1.444, 1.50, 5.0, 3.0   # rectangle: quasi-TE fundamental
    lam = 1.55
    dx = dy = 0.12
    n_ref = 0.5 * (n_core + n_clad) + 0.4 * (n_core - n_clad)
    s = BPM3DSolver(wavelength_um=lam, n_ref=n_ref, dx_um=dx, dy_um=dy, dz_um=0.5,
                    pml_um=3.0, pol="full")
    x = s.make_symmetric_grid(-11, 11, dx)
    y = s.make_symmetric_grid(-11, 11, dy)
    core = (np.abs(x)[:, None] <= wx / 2) & (np.abs(y)[None, :] <= wy / 2)
    n2 = np.where(core, n_core ** 2, n_clad ** 2)
    idx = np.sqrt(n2)

    Ex, Ey, neff = s.solve_mode_fullvec(x, y, n2, 0, 0, wx / 2.2, wy / 2.2,
                                        launch="Ex", max_iter=700)
    assert n_clad < neff < n_core, f"hybrid mode n_eff {neff} out of range"
    p_ex = _comp_power(Ex, x, y)
    p_ey = _comp_power(Ey, x, y)
    conv_in = p_ey / (p_ex + p_ey)

    z = np.arange(0, 90 + 0.25, s.dz_um)
    res = s.run_fullvec(x, y, z, lambda iz: idx, Ex, Ey, save_every=10)
    # shape fidelity: normalise the output vector field, overlap with the mode
    pout = np.sqrt(_comp_power(res.psi_out, x, y) + _comp_power(res.psi_out_minor, x, y))
    sov = s.vector_overlap_power(res.psi_out / pout, res.psi_out_minor / pout, Ex, Ey, x, y)

    # Gaussian Ex launch (not an eigenstate) breathes more
    g = s.gaussian_2d(x, y, 0, 0, wx / 2.2, wy / 2.2)
    res_g = s.run_fullvec(x, y, z, lambda iz: idx, g, np.zeros_like(g), save_every=10)
    pg = np.sqrt(_comp_power(res_g.psi_out, x, y) + _comp_power(res_g.psi_out_minor, x, y))
    sov_g = s.vector_overlap_power(res_g.psi_out / pg, res_g.psi_out_minor / pg, g, np.zeros_like(g), x, y)

    print(f"[fullvec hybrid] n_eff={neff:.5f}  minor_in={conv_in*100:.3f}%  "
          f"vec self-overlap mode={sov*100:.1f}%  gaussian={sov_g*100:.1f}%")
    assert sov > 0.96, f"hybrid mode should stay stationary, got {sov}"
    assert sov > sov_g, "hybrid mode must breathe less than a Gaussian launch"
    print("  PASS")


def _silica_n(lam):
    lam2 = lam ** 2
    b = (0.6961663, 0.4079426, 0.8974794)
    c = (0.0684043 ** 2, 0.1162414 ** 2, 9.896161 ** 2)
    n2 = 1 + sum(bi * lam2 / (lam2 - ci) for bi, ci in zip(b, c))
    return n2 ** 0.5


def test_port_modal_selectivity():
    """The per-port modal-matrix building block: isolating a port's guide and
    solving its mode must yield a mode localized at that port — overlapping a
    field guided in port 1 couples strongly to mode 1, negligibly to mode 2.
    Also sanity-check the silica dispersion ratio direction."""
    n_clad, n_core, w = 1.444, 1.49, 4.0
    dx = dy = 0.15
    n_ref = 0.5 * (n_core + n_clad) + 0.4 * (n_core - n_clad)
    s = BPM3DSolver(wavelength_um=1.55, n_ref=n_ref, dx_um=dx, dy_um=dy, dz_um=0.5, pml_um=3.0)
    x = s.make_symmetric_grid(-18, 18, dx)
    y = s.make_symmetric_grid(-12, 12, dy)
    X = 7.0
    core1 = (np.abs(x - X)[:, None] <= w / 2) & (np.abs(y)[None, :] <= w / 2)
    core2 = (np.abs(x + X)[:, None] <= w / 2) & (np.abs(y)[None, :] <= w / 2)
    n2 = np.where(core1 | core2, n_core ** 2, n_clad ** 2)
    n2c = float(n2.min())

    def port_idx(cx):                       # isolate the guide near cx
        near = np.abs(x - cx) <= w * 1.5
        return np.where(near[:, None], n2, n2c)

    m1, nf1 = s.solve_mode(x, y, port_idx(+X), +X, 0, w / 2.2, w / 2.2, max_iter=500)
    m2, nf2 = s.solve_mode(x, y, port_idx(-X), -X, 0, w / 2.2, w / 2.2, max_iter=500)
    assert n_clad < nf1 < n_core and n_clad < nf2 < n_core
    cross = s.overlap_power(m1, m2, x, y)   # disjoint guides -> ~0
    assert cross < 1e-3, f"isolated port modes should be near-orthogonal, got {cross}"

    idx = np.sqrt(n2)
    z = np.arange(0, 60 + 0.25, s.dz_um)
    res = s.run(x, y, z, lambda iz: idx, m1, save_every=10)   # launch port-1 mode
    o1 = s.overlap_power(res.psi_out, m1, x, y) / res.input_power
    o2 = s.overlap_power(res.psi_out, m2, x, y) / res.input_power
    print(f"[port modal] n_eff1={nf1:.5f} n_eff2={nf2:.5f}  cross={cross:.1e}  "
          f"coupled into port1={o1*100:.1f}% port2={o2*100:.2f}%")
    assert o1 > 0.85, f"port-1 launch should stay in port 1, got {o1}"
    assert o2 < 0.05, f"negligible power should reach port 2, got {o2}"

    # silica material dispersion: normal dispersion -> n decreases with wavelength
    assert _silica_n(1.50) > _silica_n(1.60), "silica dispersion ratio direction wrong"
    print("  PASS")


def test_gui_pipeline():
    """Exercise the exact data path run_bpm_3d uses: a (nz,nx) lateral mask +
    VerticalStack -> index_fn closure -> 3D run.  (No GDS/shapely needed.)"""
    n_clad, n_core, h = 1.444, 1.4549, 6.0
    dx, dy, dz = 0.15, 0.15, 0.5

    # common (x, z) grid like build_core_mask_grid would produce
    x = BPM3DSolver.make_symmetric_grid(-18, 18, dx)
    z = np.arange(0, 200 + 0.5 * dz, dz)
    nz, nx = len(z), len(x)

    # straight 6um core on layer (1,0)
    core_row = np.abs(x) <= h / 2
    mask = np.broadcast_to(core_row, (nz, nx)).copy()
    masks = {(1, 0): mask}

    stack = default_silica_stack(n_core, n_clad, core_thickness_um=h,
                                 gds_layer=1, gds_datatype=0)
    assert stack.mask_keys() == [(1, 0)]
    assert abs(stack.core_center_um()) < 1e-9
    assert abs(stack.index_max() - n_core) < 1e-9

    y = stack.make_y_grid(dy, margin_um=4.0)
    ny = len(y)
    bg_profile = stack.background_profile(y)
    bands = stack.patterned_bands(y)
    bg2d = np.broadcast_to(bg_profile[None, :], (nx, ny))

    def index_fn(iz):
        n_xy = bg2d.copy()
        for key, band, n_val in bands:
            lat = masks.get(key)
            row = lat[iz]
            if row.any():
                n_xy[row[:, None] & band[None, :]] = n_val
        return n_xy

    # index sanity: core box present, background elsewhere.
    # patterned_bands uses a half-open vertical band [y0, y0+h); match it here.
    n0 = index_fn(0)
    core_box = (np.abs(x)[:, None] <= h / 2) & ((y >= -h / 2) & (y < h / 2))[None, :]
    assert np.allclose(n0[core_box], n_core), "core index not painted"
    assert np.allclose(n0[~core_box], n_clad), "background index wrong"

    # two-step EIM reference index (mirrors _estimate_n_ref_3d without bpm_gui import)
    n_ref = 0.5 * (n_core + n_clad) + 0.4 * (n_core - n_clad)
    s = BPM3DSolver(wavelength_um=1.55, n_ref=n_ref, dx_um=dx, dy_um=dy, dz_um=dz, pml_um=3.0)
    field0 = s.gaussian_2d(x, y, 0, 0, h / 2.2, h / 2.2)
    res = s.run(x, y, z, index_fn, field0, save_every=10)

    win_x, win_y = np.abs(x) <= h, np.abs(y) <= h
    def conf(fld):
        sub = fld[np.ix_(win_x, win_y)]
        return float(np.trapz(np.trapz(np.abs(sub) ** 2, y[win_y], axis=1), x[win_x]))
    c0, c1 = conf(field0), conf(res.psi_out)
    p0, p1 = res.total_power[0], res.total_power[-1]
    print(f"[gui pipeline] nx={nx} ny={ny} nz={nz}  power {p0:.4f}->{p1:.4f}  "
          f"confined {c1/c0*100:.1f}%  topview{res.topview.shape} sideview{res.sideview.shape}")
    assert c1 / c0 > 0.80, "pipeline: beam not guided"
    assert abs(p1 - p0) / p0 < 0.08, "pipeline: power drift"
    assert res.topview.shape[1] == nx and res.sideview.shape[1] == ny
    print("  PASS")


if __name__ == "__main__":
    test_free_space_diffraction()
    test_straight_waveguide()
    test_pml_absorption()
    test_mode_slab_te_tm()
    test_eigenmode_launch()
    test_fullvec_separable_no_coupling()
    test_fullvec_power_conservation()
    test_fullvec_rotated_conversion()
    test_fullvec_hybrid_mode()
    test_port_modal_selectivity()
    test_gui_pipeline()
    print("\nALL TESTS PASSED")
