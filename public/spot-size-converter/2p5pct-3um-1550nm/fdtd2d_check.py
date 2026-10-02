"""
fdtd2d_check.py -- rigorous calibration of the segment-radiation loss
=====================================================================
The 3D BPM is scalar and paraxial: radiation that a z-segmented core sheds at
large angles (first grating order of a multi-um pitch leaves at ~40-60 deg) is
only approximately represented.  This script quantifies that error on a 2D
(vertical y-z) slab version of the problem -- 3 um thick 2.5% Delta core -- by
comparing, for the SAME segmented section,

  (a) a full-wave 2D FDTD (TE, Ex-Hy-Hz, absorbing layers, CW + DFT), and
  (b) the same scalar paraxial Crank-Nicolson BPM used in 3D,

each reporting the transmission into the fundamental mode of a solid output
slab, normalized to an unsegmented reference.  The ratio (b)/(a) of the excess
loss is the calibration factor for the 3D BPM segment loss.

Units: um, c = 1.
"""
import math, sys, json
import numpy as np
from numba import njit
from ssc_core import Platform

P0 = Platform()
LAM = P0.lam_um
NCL, NCO = P0.n_clad, P0.n_core
H = P0.core_h_um


def slab_mode(y, n):
    """Fundamental TE mode of a 1D slab profile n(y) (FD, Dirichlet)."""
    dy = y[1] - y[0]; k0 = 2 * math.pi / LAM
    N = len(y)
    A = (np.diag(-2 * np.ones(N)) + np.diag(np.ones(N - 1), 1) + np.diag(np.ones(N - 1), -1)) / dy ** 2
    A += np.diag((k0 * n) ** 2)
    w, v = np.linalg.eigh(A)
    m = v[:, -1]
    m = m * np.sign(m[N // 2]) / math.sqrt(np.sum(m ** 2) * dy)
    return math.sqrt(w[-1]) / k0, m


def structure(z, y, pitch, duty, L_seg, z0, solid_core=True):
    """eps(y,z): solid core, then [z0, z0+L_seg) segmented (tooth first), then solid."""
    core = np.abs(y)[None, :] <= H / 2
    ph = ((z - z0) % pitch) / pitch
    inseg = (z >= z0) & (z < z0 + L_seg)
    tooth = np.where(inseg, ph < duty, True)
    n = np.where(core & tooth[:, None], NCO, NCL)
    return n


@njit(cache=True)
def _fdtd(eps, sig, src_iz, src_prof, omega, nsteps, dt, dz, dy, dft_iz, ndft_start):
    nz, ny = eps.shape
    Ex = np.zeros((nz, ny)); Hy = np.zeros((nz, ny)); Hz = np.zeros((nz, ny))
    ca = (1 - sig * dt / (2 * eps)) / (1 + sig * dt / (2 * eps))
    cb = (dt / eps) / (1 + sig * dt / (2 * eps))
    sigm = sig / eps            # matched magnetic loss (mu = 1)
    da = (1 - sigm * dt / 2) / (1 + sigm * dt / 2)
    db_ = dt / (1 + sigm * dt / 2)
    acc_re = np.zeros(ny); acc_im = np.zeros(ny)
    for n in range(nsteps):
        t = n * dt
        # H updates
        for i in range(nz - 1):
            for j in range(ny):
                Hy[i, j] = da[i, j] * Hy[i, j] - db_[i, j] * (Ex[i + 1, j] - Ex[i, j]) / dz
        for i in range(nz):
            for j in range(ny - 1):
                Hz[i, j] = da[i, j] * Hz[i, j] + db_[i, j] * (Ex[i, j + 1] - Ex[i, j]) / dy
        # E update
        for i in range(1, nz):
            for j in range(1, ny):
                curl = (Hz[i, j] - Hz[i, j - 1]) / dy - (Hy[i, j] - Hy[i - 1, j]) / dz
                Ex[i, j] = ca[i, j] * Ex[i, j] + cb[i, j] * curl
        ramp = min(1.0, t / 60.0)
        s = ramp * math.sin(omega * t)
        for j in range(ny):
            Ex[src_iz, j] += src_prof[j] * s * dt
        if n >= ndft_start:
            c = math.cos(omega * t); sn = math.sin(omega * t)
            for j in range(ny):
                acc_re[j] += Ex[dft_iz, j] * c
                acc_im[j] += Ex[dft_iz, j] * sn
    return acc_re + 1j * acc_im


def fdtd_T(pitch, duty, L_seg, d=0.04, ylim=12.0, pml=2.5):
    z = np.arange(0, 6 + L_seg + 14 + 2 * pml, d)
    y = np.arange(-ylim, ylim + 1e-9, d)
    z0 = pml + 4.0
    n = structure(z, y, pitch, duty, L_seg, z0) if L_seg > 0 else structure(z, y, 1, 1, 0, z0)
    eps = n ** 2
    # absorbing layers
    sig = np.zeros_like(eps)
    smax = 3.0
    dzb = np.clip((pml - z) / pml, 0, 1) + np.clip((z - (z[-1] - pml)) / pml, 0, 1)
    dyb = np.clip((np.abs(y) - (ylim - pml)) / pml, 0, 1)
    sig += smax * eps * np.maximum(dzb[:, None] ** 2, dyb[None, :] ** 2)
    ncore_prof = np.where(np.abs(y) <= H / 2, NCO, NCL)
    neff, m = slab_mode(y, ncore_prof)
    omega = 2 * math.pi / LAM
    dt = 0.5 * d / math.sqrt(2)
    src_iz = int((pml + 1.0) / d)
    dft_iz = int((z0 + L_seg + 8.0) / d)
    T_total = (z[-1] * NCO) * 2.2 + 80
    nsteps = int(T_total / dt)
    ndft = int((T_total - 30 * LAM) / dt)
    Eph = _fdtd(eps, sig, src_iz, m.copy(), omega, nsteps, dt, d, d, dft_iz, ndft)
    amp = np.sum(Eph * m) * d
    return abs(amp) ** 2


def bpm2d_T(pitch, duty, L_seg, d=0.04, dz=0.05, ylim=12.0, pml=2.5):
    y = np.arange(-ylim, ylim + 1e-9, d)
    k0 = 2 * math.pi / LAM
    ncore_prof = np.where(np.abs(y) <= H / 2, NCO, NCL)
    neff, m = slab_mode(y, ncore_prof)
    b0 = k0 * (NCL + 0.004)
    N = len(y)
    D2 = (np.diag(-2 * np.ones(N)) + np.diag(np.ones(N - 1), 1) + np.diag(np.ones(N - 1), -1)) / d ** 2
    ab = 2.0 * b0 * 0.5 * np.clip((np.abs(y) - (ylim - pml)) / pml, 0, 1) ** 2
    def op(nprof):
        return D2 + np.diag(k0 ** 2 * nprof ** 2 - b0 ** 2 + 1j * ab)
    I = np.eye(N)
    a = 1j * dz / (4 * b0)          # Crank-Nicolson: (1 - a L) psi+ = (1 + a L) psi
    Ms = {}
    for key, prof in (("t", ncore_prof), ("g", np.full(N, NCL))):
        Lop = op(prof)
        Ms[key] = np.linalg.solve(I - a * Lop, I + a * Lop)
    psi = m.astype(complex)
    nst = int(round(L_seg / dz))
    for k in range(nst):
        zc = (k + 0.5) * dz
        ph = (zc % pitch) / pitch
        psi = Ms["t" if ph < duty else "g"] @ psi
    # 8 um solid output then overlap
    for k in range(int(8 / dz)):
        psi = Ms["t"] @ psi
    return abs(np.sum(psi * m) * d) ** 2


if __name__ == "__main__":
    cases = [(float(a), float(b)) for a, b in (s.split(",") for s in sys.argv[1:])] or \
            [(2.0, 0.6), (2.0, 0.4), (2.8, 0.4), (2.8, 0.25)]
    L = 60.0
    Tref_f = fdtd_T(1.0, 1.0, 0.0)
    out = []
    for P, du in cases:
        Lseg = round(L / P) * P
        Tf = fdtd_T(P, du, Lseg) / Tref_f
        Tb = bpm2d_T(P, du, Lseg)
        lf, lb = -10 * math.log10(Tf), -10 * math.log10(Tb)
        out.append(dict(pitch=P, duty=du, L=Lseg, fdtd_dB=lf, bpm_dB=lb))
        print(f"P={P:.2f} duty={du:.2f} L={Lseg:.1f}um  FDTD {lf:.4f} dB   paraxial BPM {lb:.4f} dB  "
              f"(per-period: {lf/(Lseg/P):.5f} vs {lb/(Lseg/P):.5f})", flush=True)
    json.dump(out, open("fdtd2d_calibration.json", "w"), indent=2)
