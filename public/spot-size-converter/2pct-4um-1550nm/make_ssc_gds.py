"""Build a complete SSC device GDS from the 2% 4x4um @1550nm page result.

Page design (mode-overlap, 0.114 dB to SMF-28 @1550nm):
  chip 4x4um solid core (2% Delta GeO2-silica)
  mode-expanding segmented taper: pitch 4.0um, 120 segments,
      duty 0.95 -> 0.20 (cosine),  width 4.0 -> 7.0um (cosine)
  facet (SMF butt-coupling): 7.0um wide, 0.20 duty  (expanded mode ~11um)

Layout (centred on y = 0, propagation along +z):
  layer (1,0)  core  : solid input bus + segmented taper (teeth) + facet
  layer (10,0) facet : dicing / SMF butt-coupling reference line at the chip edge
  layer (63,0) text  : design label
"""
import numpy as np
import gdstk
from seg_ssc import SegSSCGene

# --- exact page design ---
gene = SegSSCGene(pitch_um=4.0, n_seg=120, duty_start=0.95, duty_end=0.20,
                  w_start=4.0, w_end=7.0, duty_profile="cos", width_profile="cos",
                  leadins=[(50.0, 4.0)])

CORE_LAYER = (1, 0)
FACET_LAYER = (10, 0)
TEXT_LAYER = (63, 0)
BUS_IN_UM = 150.0          # solid 4um input bus routed onto the chip (z < 0)
W_CHIP = gene.w_start      # 4.0 um

lib = gdstk.Library(unit=1e-6, precision=1e-9)
cell = lib.new_cell("SSC_2PCT_4UM_1550")

# (1) solid input chip bus (z from -BUS_IN to 0), 4 um wide
cell.add(gdstk.rectangle((-BUS_IN_UM, -W_CHIP / 2), (0.0, W_CHIP / 2),
                         layer=CORE_LAYER[0], datatype=CORE_LAYER[1]))

# (2) segmented mode-expanding taper: one rectangle per solid tooth
for (z0, z1, w) in gene.segments():
    hw = 0.5 * w
    cell.add(gdstk.rectangle((z0, -hw), (z1, hw),
                             layer=CORE_LAYER[0], datatype=CORE_LAYER[1]))

L = gene.total_length      # facet z-position (chip edge)

# (3) facet / dicing reference line at the chip edge (SMF butt-couples here)
FW = 0.5 * gene.w_end + 6.0
cell.add(gdstk.rectangle((L - 0.2, -FW), (L + 0.2, FW),
                         layer=FACET_LAYER[0], datatype=FACET_LAYER[1]))

# (4) design label
label = ("2% 4x4um SSC @1550nm  |  chip 4um -> facet 7.0um/duty0.20  |  "
         "L={:.0f}um  |  SMF-28 coupling 0.114 dB (mode-overlap)").format(L)
cell.add(gdstk.Label(label, (L * 0.5, FW + 3.0), layer=TEXT_LAYER[0],
                     texttype=TEXT_LAYER[1], magnification=4.0))

lib.write_gds("ssc_device.gds")

# --- summary ---
segs = gene.segments()
print(f"cell: {cell.name}")
print(f"total taper length : {L:.1f} um  (+ {BUS_IN_UM:.0f} um input bus)")
print(f"segments (incl. lead-in): {len(segs)}")
print(f"pitch / n_seg      : {gene.pitch_um} um / {gene.n_seg}")
print(f"duty  start->end   : {gene.duty_start} -> {gene.duty_end} ({gene.duty_profile})")
print(f"width start->end   : {gene.w_start} -> {gene.w_end} um ({gene.width_profile})")
bb = cell.bounding_box()
print(f"bbox (um)          : {tuple(round(v,2) for v in bb[0])} -> {tuple(round(v,2) for v in bb[1])}")
print(f"polygons           : {len(cell.polygons)}   layers: {sorted({(p.layer,p.datatype) for p in cell.polygons})}")
print("wrote ssc_device.gds")

# --- top-view preview PNG ---
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

fig, ax = plt.subplots(2, 1, figsize=(11, 4.2), gridspec_kw={"height_ratios": [2, 1]})
# full device
ax[0].add_patch(Rectangle((-BUS_IN_UM, -W_CHIP/2), BUS_IN_UM, W_CHIP, color="#3b5bdb"))
for (z0, z1, w) in segs:
    ax[0].add_patch(Rectangle((z0, -w/2), z1 - z0, w, color="#3b5bdb"))
ax[0].axvline(L, color="#e03131", lw=1.2, ls="--")
ax[0].text(L, 7.5, "facet / SMF", color="#e03131", ha="center", fontsize=8)
ax[0].set_xlim(-BUS_IN_UM, L + 20); ax[0].set_ylim(-8, 10)
ax[0].set_xlabel("z (propagation) [um]"); ax[0].set_ylabel("x [um]")
ax[0].set_title("2% 4×4µm SSC device @1550 nm — full layout (core layer)")
ax[0].set_aspect("auto")
# zoom on the facet (last ~60um) to show the dashed teeth
z0z = L - 60
for (a, b, w) in segs:
    if b > z0z:
        ax[1].add_patch(Rectangle((a, -w/2), b - a, w, color="#3b5bdb"))
ax[1].axvline(L, color="#e03131", lw=1.2, ls="--")
ax[1].set_xlim(z0z, L + 5); ax[1].set_ylim(-5, 5)
ax[1].set_xlabel("z [um]"); ax[1].set_ylabel("x [um]")
ax[1].set_title("Facet-end zoom (low-duty segments)")
ax[1].set_aspect("equal")
fig.tight_layout()
fig.savefig("ssc_device_layout.png", dpi=130)
print("wrote ssc_device_layout.png")
