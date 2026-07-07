"""Robust 2-D finite-difference scalar/semi-vectorial mode solver for the
Si3N4 (in SiO2) inverse-taper SSC.  Uses a sparse shift-invert eigensolve,
which does not suffer the imaginary-distance BPM's spurious-mode collapse for
high-contrast waveguides near cutoff.

Scalar Helmholtz:  (d2/dx2 + d2/dy2) E + k0^2 n(x,y)^2 E = beta^2 E
Find the eigenpair with the largest beta (n_eff).  Dirichlet edges (mode is
localised; window is padded well beyond the mode).
"""
import numpy as np, math
import scipy.sparse as sp
import scipy.sparse.linalg as spla

LAM = 1.55
def n_si3n4(l):
    l2=l*l; return math.sqrt(1+3.0249*l2/(l2-0.1353406**2)+40314.0*l2/(l2-1239.842**2))
def n_sio2(l):
    l2=l*l; return math.sqrt(1+0.6961663*l2/(l2-0.0684043**2)+0.4079426*l2/(l2-0.1162414**2)+0.8974794*l2/(l2-9.896161**2))
N_SIN = n_si3n4(LAM)
N_OX  = n_sio2(LAM)


def build_index(x, y, w_um, h_um=0.30, edge=0.015):
    """Si3N4 core (w x h) centred at origin, buried in SiO2 (soft tanh edges)."""
    hw, ht = 0.5*w_um, 0.5*h_um
    lat = 0.5*(np.tanh((x+hw)/edge) - np.tanh((x-hw)/edge))
    ver = 0.5*(np.tanh((y+ht)/edge) - np.tanh((y-ht)/edge))
    n = N_OX + (N_SIN - N_OX)*np.outer(lat, ver)   # (nx, ny)
    return n


def solve_fd(x, y, n, k0, num=1):
    """Largest-neff scalar mode(s). Returns (neff_list, modes list of (nx,ny))."""
    nx, ny = len(x), len(y)
    dx = x[1]-x[0]; dy = y[1]-y[0]
    N = nx*ny
    # 5-point Laplacian (Dirichlet), row-major index k = ix*ny + iy
    ex = np.ones(nx)/dx**2
    ey = np.ones(ny)/dy**2
    Lx = sp.diags([ex[:-1], -2*ex, ex[:-1]], [-1,0,1], shape=(nx,nx))
    Ly = sp.diags([ey[:-1], -2*ey, ey[:-1]], [-1,0,1], shape=(ny,ny))
    Lap = sp.kron(Lx, sp.identity(ny)) + sp.kron(sp.identity(nx), Ly)
    diagn = sp.diags((k0**2 * (n.reshape(N)**2)))
    A = (Lap + diagn).tocsc()
    sigma = (k0*N_SIN)**2 * 0.999   # shift near the max possible beta^2
    vals, vecs = spla.eigsh(A, k=num, sigma=sigma, which='LM')
    order = np.argsort(vals)[::-1]
    neffs=[]; modes=[]
    for idx in order:
        beta2 = vals[idx]
        neffs.append(math.sqrt(beta2)/k0 if beta2>0 else 0.0)
        modes.append(vecs[:,idx].reshape(nx,ny))
    return neffs, modes


def d4sigma(f, coord, axis):
    I=np.abs(f)**2; p=I.sum(axis=1-axis); p=p/p.sum(); m=(coord*p).sum()
    return 4*math.sqrt(((coord-m)**2*p).sum())


def gauss2d(x, y, w0):
    X,Y=np.meshgrid(x,y,indexing='ij')
    g=np.exp(-(X**2+Y**2)/w0**2); return g


def overlap(a, b, x, y):
    dx=x[1]-x[0]; dy=y[1]-y[0]
    num=abs(np.sum(a*np.conj(b)))**2
    da=np.sum(np.abs(a)**2); db=np.sum(np.abs(b)**2)
    return float(num/(da*db))


if __name__ == "__main__":
    k0 = 2*math.pi/LAM
    MFD = 4.1; w0 = MFD/2
    dx=0.035; dy=0.025
    x=np.arange(-8,8+1e-9,dx); y=np.arange(-6,6+1e-9,dy)
    uhna = gauss2d(x,y,w0)
    print(f"n_SiN={N_SIN:.5f} n_SiO2={N_OX:.5f}  grid {len(x)}x{len(y)}")
    print("w_tip  neff       D4x   D4y   eta(%)  loss(dB)")
    for w in [0.14,0.16,0.18,0.20,0.22,0.24,0.26,0.30,0.40,1.50]:
        n=build_index(x,y,w)
        neffs,modes=solve_fd(x,y,n,k0,num=1)
        m=modes[0]; ne=neffs[0]
        eta=overlap(m,uhna,x,y)
        gd=(d4sigma(m,x,0), d4sigma(m,y,1))
        print(f"{w:5.2f}  {ne:.6f}  {gd[0]:4.2f}  {gd[1]:4.2f}  {eta*100:6.2f}  {-10*math.log10(max(eta,1e-9)):.4f}", flush=True)
