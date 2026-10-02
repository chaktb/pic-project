"""
ssc_core.py -- 2.5% Delta, 3 x 3 um silica core  <->  SMF-28 @ 1550 nm
=====================================================================
Self-contained numerics for the segmented spot-size converter (SSC):

  * Platform           : Sellmeier silica cladding, Delta -> n_core
  * Grid / index model : sub-pixel (area-averaged) rectangular core,
                         segmented section = duty-averaged permittivity
                         (effective-medium) or true tooth/gap along z
  * mode_solve()       : scalar finite-difference eigenmode solver
  * BPM (ADI-CN)       : scalar paraxial 3D BPM, Peaceman-Rachford ADI,
                         Crank-Nicolson, numba-batched Thomas solves.
                         It uses the SAME 5-point Laplacian and n^2 map as
                         the mode solver, so an eigenmode propagates without
                         the spurious "mode-beating" of mismatched operators.
  * SSCDesign          : geometry gene (lead-in, pitch, duty/width ramps)
                         -> z-sampled (width, duty) or real tooth list
                         with DRC (tooth, gap, width >= 0.5 um).

Delta definition: Delta = (n1^2 - n2^2) / (2 n1^2)  ->  n1 = n2 / sqrt(1 - 2 Delta)
Units: um.
"""
from __future__ import annotations
import math
from dataclasses import dataclass, field
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from numba import njit, prange

# ---------------------------------------------------------------------------
# 1. Platform
# ---------------------------------------------------------------------------

def n_fused_silica(lam_um: float) -> float:
    """Malitson (1965) Sellmeier for fused silica."""
    l2 = lam_um ** 2
    n2 = 1 + 0.6961663 * l2 / (l2 - 0.0684043 ** 2) \
           + 0.4079426 * l2 / (l2 - 0.1162414 ** 2) \
           + 0.8974794 * l2 / (l2 - 9.896161 ** 2)
    return math.sqrt(n2)


@dataclass
class Platform:
    lam_um: float = 1.550
    delta: float = 0.025          # 2.5 %
    core_w_um: float = 3.0        # chip waveguide width
    core_h_um: float = 3.0        # core thickness (fixed by deposition)
    smf_mfd_um: float = 10.4      # SMF-28 MFD @ 1550 nm
    drc_min_um: float = 0.5       # min linewidth / gap
    n_clad: float = field(init=False)
    n_core: float = field(init=False)

    def __post_init__(self):
        self.n_clad = n_fused_silica(self.lam_um)
        self.n_core = self.n_clad / math.sqrt(1.0 - 2.0 * self.delta)

    @property
    def k0(self):
        return 2 * math.pi / self.lam_um


# ---------------------------------------------------------------------------
# 2. Grid and index cross-sections
# ---------------------------------------------------------------------------

class Grid:
    def __init__(self, half_x=20.0, half_y=18.0, dx=0.25, dy=None):
        dy = dx if dy is None else dy
        self.dx, self.dy = dx, dy
        self.x = np.arange(-half_x, half_x + 1e-9, dx)
        self.y = np.arange(-half_y, half_y + 1e-9, dy)
        self.nx, self.ny = len(self.x), len(self.y)
        self.X, self.Y = np.meshgrid(self.x, self.y, indexing="ij")

    def coverage(self, coord, d, half):
        """Fraction of each cell [c-d/2, c+d/2] inside [-half, half]."""
        lo = np.maximum(coord - d / 2, -half)
        hi = np.minimum(coord + d / 2, half)
        return np.clip(hi - lo, 0, None) / d

    def n2_rect(self, plat: Platform, width, duty=1.0, height=None):
        """Permittivity map: rectangle (width x height) whose permittivity
        contrast is scaled by `duty` (effective medium of a z-segmented core)."""
        h = plat.core_h_um if height is None else height
        cx = self.coverage(self.x, self.dx, width / 2)
        cy = self.coverage(self.y, self.dy, h / 2)
        f = np.outer(cx, cy) * duty
        return plat.n_clad ** 2 + f * (plat.n_core ** 2 - plat.n_clad ** 2)

    def gaussian(self, w0):
        g = np.exp(-(self.X ** 2 + self.Y ** 2) / w0 ** 2).astype(complex)
        return self.normalize(g)

    def normalize(self, E):
        p = np.sum(np.abs(E) ** 2) * self.dx * self.dy
        return E / math.sqrt(p)

    def power(self, E):
        return float(np.sum(np.abs(E) ** 2) * self.dx * self.dy)

    def overlap(self, E1, E2):
        """Power coupling |<E1|E2>|^2 / (<E1|E1><E2|E2>)."""
        num = abs(np.sum(np.conj(E1) * E2)) ** 2
        return float(num / (np.sum(abs(E1) ** 2) * np.sum(abs(E2) ** 2)))

    def overlap_raw(self, Ein, Emode_normalized):
        """Power from (un-normalized) field Ein into a unit-power mode."""
        return float(abs(np.sum(np.conj(Emode_normalized) * Ein) * self.dx * self.dy) ** 2)


def laplacian(grid: Grid):
    def d2(n, h):
        e = np.ones(n)
        return sp.diags([e[:-1], -2 * e, e[:-1]], [-1, 0, 1]) / h ** 2
    Ix, Iy = sp.identity(grid.nx), sp.identity(grid.ny)
    return (sp.kron(d2(grid.nx, grid.dx), Iy) + sp.kron(Ix, d2(grid.ny, grid.dy))).tocsc()


_LAP_CACHE = {}


def mode_solve(grid: Grid, plat: Platform, n2):
    """Scalar FD fundamental mode.  Returns (n_eff, E normalized to unit power)."""
    key = (grid.nx, grid.ny, grid.dx, grid.dy)
    if key not in _LAP_CACHE:
        _LAP_CACHE[key] = laplacian(grid)
    L = _LAP_CACHE[key]
    k0 = plat.k0
    A = L + sp.diags((k0 ** 2) * n2.ravel())
    sigma = (k0 * plat.n_core) ** 2
    vals, vecs = spla.eigsh(A, k=1, sigma=sigma, which="LM")
    beta = math.sqrt(vals[0])
    E = vecs[:, 0].reshape(grid.nx, grid.ny)
    E = E * np.sign(E[grid.nx // 2, grid.ny // 2] or 1.0)
    return beta / k0, grid.normalize(E.astype(complex))


def d4sigma(grid: Grid, E):
    I = np.abs(E) ** 2
    P = I.sum()
    mx = (I * grid.X).sum() / P
    my = (I * grid.Y).sum() / P
    sx = math.sqrt((I * (grid.X - mx) ** 2).sum() / P)
    sy = math.sqrt((I * (grid.Y - my) ** 2).sum() / P)
    return 4 * sx, 4 * sy


# ---------------------------------------------------------------------------
# 3. ADI Crank-Nicolson scalar BPM
# ---------------------------------------------------------------------------

@njit(cache=True, parallel=True)
def _adi_x(psi, n2k2, b0sq, a, idx2, idy2, absorb):
    """One half-step: implicit in x, explicit in y.
    (1 - a Lx) out = (1 + a Ly) psi ;   L = D2 + 0.5 (k^2 n^2 - b0^2) + i*absorb/2
    """
    nx, ny = psi.shape
    out = np.empty_like(psi)
    for j in prange(ny):
        rhs = np.empty(nx, dtype=np.complex128)
        for i in range(nx):
            v = 0.5 * (n2k2[i, j] - b0sq) + 0.5j * absorb[i, j]
            lap = -2.0 * psi[i, j] * idy2
            if j > 0:
                lap += psi[i, j - 1] * idy2
            if j < ny - 1:
                lap += psi[i, j + 1] * idy2
            rhs[i] = psi[i, j] + a * (lap + v * psi[i, j])
        # tridiagonal: sub = -a*idx2, diag = 1 - a*(-2 idx2 + v), sup = -a*idx2
        cp = np.empty(nx, dtype=np.complex128)
        dp = np.empty(nx, dtype=np.complex128)
        off = -a * idx2
        for i in range(nx):
            v = 0.5 * (n2k2[i, j] - b0sq) + 0.5j * absorb[i, j]
            diag = 1.0 - a * (-2.0 * idx2 + v)
            if i == 0:
                cp[i] = off / diag
                dp[i] = rhs[i] / diag
            else:
                m = diag - off * cp[i - 1]
                cp[i] = off / m
                dp[i] = (rhs[i] - off * dp[i - 1]) / m
        out[nx - 1, j] = dp[nx - 1]
        for i in range(nx - 2, -1, -1):
            out[i, j] = dp[i] - cp[i] * out[i + 1, j]
    return out


@njit(cache=True, parallel=True)
def _adi_y(psi, n2k2, b0sq, a, idx2, idy2, absorb):
    """Second half-step: implicit in y, explicit in x."""
    nx, ny = psi.shape
    out = np.empty_like(psi)
    for i in prange(nx):
        rhs = np.empty(ny, dtype=np.complex128)
        for j in range(ny):
            v = 0.5 * (n2k2[i, j] - b0sq) + 0.5j * absorb[i, j]
            lap = -2.0 * psi[i, j] * idx2
            if i > 0:
                lap += psi[i - 1, j] * idx2
            if i < nx - 1:
                lap += psi[i + 1, j] * idx2
            rhs[j] = psi[i, j] + a * (lap + v * psi[i, j])
        cp = np.empty(ny, dtype=np.complex128)
        dp = np.empty(ny, dtype=np.complex128)
        off = -a * idy2
        for j in range(ny):
            v = 0.5 * (n2k2[i, j] - b0sq) + 0.5j * absorb[i, j]
            diag = 1.0 - a * (-2.0 * idy2 + v)
            if j == 0:
                cp[j] = off / diag
                dp[j] = rhs[j] / diag
            else:
                m = diag - off * cp[j - 1]
                cp[j] = off / m
                dp[j] = (rhs[j] - off * dp[j - 1]) / m
        out[i, ny - 1] = dp[ny - 1]
        for j in range(ny - 2, -1, -1):
            out[i, j] = dp[j] - cp[j] * out[i, j + 1]
    return out


class BPM:
    """psi_z = (i / 2 b0) [ Dxx + Dyy + k0^2 n^2 - b0^2 ] psi   (paraxial, scalar)."""

    def __init__(self, grid: Grid, plat: Platform, n_ref=None, pml_um=3.0, sigma_max=None):
        self.g, self.p = grid, plat
        self.k0 = plat.k0
        n_ref = plat.n_clad + 0.004 if n_ref is None else n_ref
        self.b0 = self.k0 * n_ref
        # absorber: smooth quadratic ramp of an imaginary k^2 n^2 term
        sigma_max = 2.0 * self.b0 * 0.5 if sigma_max is None else sigma_max  # [1/um] -> amplitude decay rate sigma/(2 b0)
        def ramp(c, half):
            d = np.clip((np.abs(c) - (half - pml_um)) / pml_um, 0, 1)
            return d ** 2
        rx = ramp(grid.x, grid.x[-1])
        ry = ramp(grid.y, grid.y[-1])
        self.absorb = sigma_max * np.maximum(rx[:, None], ry[None, :])
        self.idx2, self.idy2 = 1 / grid.dx ** 2, 1 / grid.dy ** 2

    def step(self, psi, n2, dz):
        a = 1j * dz / (4 * self.b0)
        n2k2 = n2 * self.k0 ** 2
        b0sq = self.b0 ** 2
        psi = _adi_x(psi, n2k2, b0sq, a, self.idx2, self.idy2, self.absorb)
        psi = _adi_y(psi, n2k2, b0sq, a, self.idx2, self.idy2, self.absorb)
        return psi

    def run(self, psi, slices, record_every=None):
        """slices: iterable of (n2_map, dz, n_steps).  Returns final psi and an
        optional record of |psi|^2 on the y=0 (x-z) and x=0 (y-z) cuts."""
        recx, recy, zs, z = [], [], [], 0.0
        jy0 = self.g.ny // 2
        ix0 = self.g.nx // 2
        k = 0
        for n2, dz, nst in slices:
            for _ in range(nst):
                psi = self.step(psi, n2, dz)
                z += dz
                k += 1
                if record_every and k % record_every == 0:
                    recx.append(np.abs(psi[:, jy0]) ** 2)
                    recy.append(np.abs(psi[ix0, :]) ** 2)
                    zs.append(z)
        rec = (np.array(zs), np.array(recx), np.array(recy)) if record_every else None
        return psi, rec


# ---------------------------------------------------------------------------
# 4. SSC geometry gene
# ---------------------------------------------------------------------------

@dataclass
class SSCDesign:
    """Chip WG (W0 solid) -> solid lead -> segmented taper -> SMF facet.

    Segmented section: pitch P constant, N = round(L/P) periods.  In period k
    (s = (k+0.5)/N) the duty and width follow
        duty(s)  = d0 + (d1 - d0) * s**p_d
        width(s) = W0 + (W1 - W0) * s**p_w
    with d0 = 1 - g_min/P (first gap = DRC minimum) and d1 the facet duty.
    Each period = tooth (duty*P) followed by gap ((1-duty)*P); the facet is
    cut at the end of the last TOOTH (the final gap is omitted), so the fiber
    sees the last core segment.
    """
    L_um: float = 800.0         # segmented taper length
    pitch_um: float = 3.0
    d_facet: float = 0.25
    w_facet_um: float = 6.0
    p_d: float = 1.0
    p_w: float = 1.0
    d_start: float | None = None   # default: 1 - drc/P
    lead_um: float = 20.0          # solid 3 um lead before the taper
    W0_um: float = 3.0
    drc_um: float = 0.5

    def __post_init__(self):
        P = self.pitch_um
        if self.d_start is None:
            self.d_start = 1.0 - self.drc_um / P
        # enforce DRC on the facet tooth
        self.d_facet = float(np.clip(self.d_facet, self.drc_um / P, self.d_start))
        self.w_facet_um = max(self.w_facet_um, self.drc_um)

    @property
    def n_per(self):
        return max(2, int(round(self.L_um / self.pitch_um)))

    def periods(self):
        N = self.n_per
        s = (np.arange(N) + 0.5) / N
        duty = self.d_start + (self.d_facet - self.d_start) * s ** self.p_d
        width = self.W0_um + (self.w_facet_um - self.W0_um) * s ** self.p_w
        return duty, width

    def teeth(self, grid_snap=0.001):
        """List of (z0, z1, w0, w1) trapezoidal core segments (incl. lead)."""
        duty, width = self.periods()
        P, N = self.pitch_um, self.n_per
        segs = [(0.0, self.lead_um, self.W0_um, self.W0_um)]
        z = self.lead_um
        # width at tooth edges: linear interpolation of the period-center law
        sc = (np.arange(N) + 0.5) / N
        def w_at(s):
            return self.W0_um + (self.w_facet_um - self.W0_um) * np.clip(s, 0, 1) ** self.p_w
        for k in range(N):
            t = duty[k] * P
            # snap tooth length so DRC margins hold on the manufacturing grid
            t = round(t / grid_snap) * grid_snap
            g = round((P - t) / grid_snap) * grid_snap
            s0, s1 = k / N, (k + duty[k]) / N
            segs.append((z, z + t, float(w_at(s0)), float(w_at(s1))))
            z += t
            if k < N - 1:
                z += g
        return segs

    def drc(self):
        segs = self.teeth()[1:]
        tooth = min(s[1] - s[0] for s in segs)
        gaps = [segs[i + 1][0] - segs[i][1] for i in range(len(segs) - 1)]
        gap0 = segs[0][0] - self.lead_um
        widths = min(min(s[2], s[3]) for s in segs)
        return dict(min_tooth=tooth, min_gap=min(gaps + ([gap0] if gap0 > 0 else [])),
                    min_width=widths, length=segs[-1][1])

    # effective-medium z-samples for fast BPM / EME
    def em_track(self, dz):
        """Continuous effective-medium track: list of (z_len, width, duty)."""
        out = [(self.lead_um, self.W0_um, 1.0)]
        N = self.n_per
        Ltot = N * self.pitch_um
        nst = max(1, int(round(Ltot / dz)))
        s = (np.arange(nst) + 0.5) / nst
        duty = self.d_start + (self.d_facet - self.d_start) * s ** self.p_d
        width = self.W0_um + (self.w_facet_um - self.W0_um) * s ** self.p_w
        for d, w in zip(duty, width):
            out.append((Ltot / nst, float(w), float(d)))
        return out

    def to_dict(self):
        return dict(L_um=self.L_um, pitch_um=self.pitch_um, d_start=self.d_start,
                    d_facet=self.d_facet, w_facet_um=self.w_facet_um, p_d=self.p_d,
                    p_w=self.p_w, lead_um=self.lead_um, W0_um=self.W0_um, n_periods=self.n_per)


def db(eta):
    return -10 * math.log10(max(eta, 1e-15))


# ---------------------------------------------------------------------------
# 5. True tooth/gap (segmented) slicing for 3D BPM
# ---------------------------------------------------------------------------

def segmented_slices(des: SSCDesign, grid: Grid, plat: Platform, dz_max=0.25, n_sub=2):
    """Slices for BPM.run() following the real GDS teeth.  Each trapezoidal
    tooth is split into n_sub constant-width pieces; teeth/gaps are stepped
    with dz <= dz_max.  The gap index map (pure cladding) is shared."""
    n2_clad = np.full((grid.nx, grid.ny), plat.n_clad ** 2)
    out = []
    segs = des.teeth()
    zprev = 0.0
    for (z0, z1, w0, w1) in segs:
        if z0 > zprev + 1e-9:                       # gap
            gl = z0 - zprev
            n = max(1, int(math.ceil(gl / dz_max)))
            out.append((n2_clad, gl / n, n))
        sub = 1 if abs(w1 - w0) < 1e-6 else n_sub
        for s in range(sub):
            za = z0 + (z1 - z0) * s / sub
            zb = z0 + (z1 - z0) * (s + 1) / sub
            w = w0 + (w1 - w0) * (s + 0.5) / sub
            n2 = grid.n2_rect(plat, w, 1.0)
            n = max(1, int(math.ceil((zb - za) / dz_max)))
            out.append((n2, (zb - za) / n, n))
        zprev = z1
    return out
