"""
Silica 2% Delta 20 ns Spiral Delay Line - GDSFactory generator
=============================================================
Platform          : Ge-doped SiO2 core (2% index contrast) buried channel, PLC
Design wavelength : 1550 nm
Target            : 20 ns optical group delay

Group index (material, exact)
    Fused silica Sellmeier (Malitson 1965):
        n^2 = 1 + 0.6961663 L^2/(L^2-0.0684043^2)
                + 0.4079426 L^2/(L^2-0.1162414^2)
                + 0.8974794 L^2/(L^2-9.896161^2)        [L in um]
    n_g = n - L dn/dL  ->  n(1550) = 1.4440,  n_g(1550) = 1.4626

For a 2% Delta Ge-doped core the phase index rises by ~2% over pure
silica, but the material *group* index is dominated by silica
dispersion, so n_g ~= 1.4626 is used as the design value. The true
waveguide group index depends on modal confinement of the ~6x6 um core
and shifts this slightly; scale N_LOOPS if a measured n_g is available.

Group delay -> physical length
    L = c * tau / n_g = (2.99792458e8 * 20e-9) / 1.4626 = 4.099 m

Low-contrast silica needs large bend radii (>= 1.5 mm) to avoid
radiation loss, so the ~4.10 m path folds into a ~11.8 x 11.8 mm die.
Continuous in-and-out double (Archimedean) spiral, ports at the center,
no crossing.
"""
import gdsfactory as gf

gf.gpdk.PDK.activate()

# --- design point ------------------------------------------------------
LAMBDA_NM = 1550.0
C_LIGHT   = 2.99792458e8   # m/s
TAU       = 20e-9          # s   (20 ns)

# exact material group index from fused-silica Sellmeier (Malitson 1965)
def _n_silica(lam_um):
    l2 = lam_um * lam_um
    return (1 + 0.6961663*l2/(l2-0.0684043**2)
              + 0.4079426*l2/(l2-0.1162414**2)
              + 0.8974794*l2/(l2-9.896161**2))**0.5
def _group_index(nfunc, lam_um, dl=1e-4):
    n0 = nfunc(lam_um)
    dndl = (nfunc(lam_um+dl) - nfunc(lam_um-dl)) / (2*dl)
    return n0 - lam_um*dndl

N_PHASE = _n_silica(LAMBDA_NM/1000)               # 1.4440
N_GROUP = _group_index(_n_silica, LAMBDA_NM/1000) # 1.4626
L_TARGET_M = C_LIGHT * TAU / N_GROUP

# --- geometry ----------------------------------------------------------
WG_WIDTH   = 6.0      # um   core width (2% Delta, ~6x6 um single mode)
MIN_RADIUS = 1500.0   # um   minimum bend radius (low-loss for 2% silica)
SEPARATION = 25.0     # um   pitch between adjacent turns
N_LOOPS    = 88       # solved so realised length ~= 20 ns target

xs = gf.cross_section.strip(width=WG_WIDTH, radius=MIN_RADIUS)
c = gf.components.spiral_double(
    min_bend_radius=MIN_RADIUS, separation=SEPARATION,
    number_of_loops=N_LOOPS, npoints=8000, cross_section=xs,
)

# --- report ------------------------------------------------------------
L_real_um = c.info["length"]
bb = c.bbox(); w, h = bb.right-bb.left, bb.top-bb.bottom
print(f"Wavelength      : {LAMBDA_NM:.0f} nm")
print(f"n (phase)       : {N_PHASE:.4f}")
print(f"n_g (group)     : {N_GROUP:.4f}")
print(f"Target delay    : {TAU*1e9:.1f} ns")
print(f"Target length   : {L_TARGET_M*1e3:.1f} mm")
print(f"Realised length : {L_real_um/1e6:.4f} m")
print(f"Realised delay  : {L_real_um*1e-6*N_GROUP/C_LIGHT*1e9:.2f} ns")
print(f"Die footprint   : {w/1000:.2f} x {h/1000:.2f} mm")

c.write_gds("silica_2pct_20ns_delay.gds")
print("Written: silica_2pct_20ns_delay.gds")
