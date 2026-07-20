"""
Si3N4 20 ns Spiral Delay Line - GDSFactory generator
====================================================
Platform : LPCVD stoichiometric Si3N4 strip waveguide (SiO2 clad)
Target   : 20 ns optical group delay @ 1550 nm

Group delay -> physical length
    L = c * tau / n_g = (2.99792458e8 * 20e-9) / 2.0 = 2.998 m

n_g = 2.0 : TE00 waveguide group index of a ~1.0 um wide Si3N4 core at
            1550 nm. Si3N4 tolerates tight bends, so a 100 um minimum
            radius keeps radiation loss negligible while the 3 m path
            folds into a ~7.6 x 7.6 mm die.

The path is a continuous in-and-out double (Archimedean) spiral: input
and output ports sit on the same edge with no crossing. The loop count
is chosen so the realised path length lands on the 20 ns target.
"""
import gdsfactory as gf

gf.gpdk.PDK.activate()  # generic PDK required for the spiral cells

# --- physical target ---------------------------------------------------
C_LIGHT = 2.99792458e8   # m/s
TAU     = 20e-9          # s   (20 ns)
N_GROUP = 2.0            # Si3N4 TE00 waveguide group index @1550 nm
L_TARGET_M  = C_LIGHT * TAU / N_GROUP
L_TARGET_UM = L_TARGET_M * 1e6

# --- geometry ----------------------------------------------------------
WG_WIDTH   = 1.0     # um   core width
MIN_RADIUS = 100.0   # um   minimum bend radius (low-loss for Si3N4)
SEPARATION = 15.0    # um   pitch between adjacent turns
N_LOOPS    = 123     # solved so realised length ~= 20 ns target

xs = gf.cross_section.strip(width=WG_WIDTH, radius=MIN_RADIUS)

c = gf.components.spiral_double(
    min_bend_radius=MIN_RADIUS,
    separation=SEPARATION,
    number_of_loops=N_LOOPS,
    npoints=8000,
    cross_section=xs,
)


# --- report ------------------------------------------------------------
L_real_um = c.info["length"]
bb = c.bbox()
w, h = bb.right - bb.left, bb.top - bb.bottom
print(f"Target delay    : {TAU*1e9:.1f} ns  (n_g = {N_GROUP})")
print(f"Target length   : {L_TARGET_M*1e3:.1f} mm")
print(f"Realised length : {L_real_um/1e6:.4f} m  "
      f"(delay = {L_real_um*1e-6*N_GROUP/C_LIGHT*1e9:.2f} ns)")
print(f"Die footprint   : {w/1000:.2f} x {h/1000:.2f} mm")

c.write_gds("si3n4_20ns_delay.gds")
print("Written: si3n4_20ns_delay.gds")
