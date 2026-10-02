"""
make_gds.py -- GDS for the GA-optimized 2.5% Delta 3x3 um segmented SSC @ 1550 nm
=================================================================================
Cells
  SSC_2p5_3um_1550        single SSC: solid 3 um lead -> segmented duty taper -> facet (+x end)
  SSC_B2B_TEST            SSC - 3 mm straight 3 um WG - SSC (back-to-back, for cut-back / IL test)
  REF_STRAIGHT            bare 3 um straight WG of the same total length (no SSC)
  SSC_ALT_INV_TAPER       alternative: solid inverse taper 3 -> 0.55 um, 2 mm (no gaps)
  TOP                     the four above, stacked with 250 um pitch

Layers
  (1, 0)   WG core (2.5% Delta GeO2-SiO2, 3 um thick)
  (2, 0)   FACET     dicing / polishing line marker at the SSC facet
  (10, 0)  TEXT      labels (non-fabrication)

Geometry is snapped to a 1 nm grid; the DRC report checks tooth length,
gap length and linewidth >= 0.5 um.
"""
import json, os, sys
import numpy as np
import gdstk
from ssc_core import SSCDesign

HERE = os.path.dirname(os.path.abspath(__file__))
L_WG, L_TEXT, L_FACET = (1, 0), (10, 0), (2, 0)
GRID = 0.001


def load_design(path=None):
    for p in ([path] if path else []) + ["final_design.json", "ga_stage3_refine.json", "ga_stage2_seg.json"]:
        fp = os.path.join(HERE, p)
        if p and os.path.exists(fp):
            d = json.load(open(fp))
            dd = d["design"]
            return SSCDesign(L_um=dd["n_periods"] * dd["pitch_um"], pitch_um=dd["pitch_um"],
                             d_facet=dd["d_facet"], w_facet_um=dd["w_facet_um"], p_d=dd["p_d"],
                             p_w=dd["p_w"], d_start=dd["d_start"], lead_um=dd["lead_um"],
                             W0_um=dd["W0_um"]), p
    raise FileNotFoundError("no design json")


def snap(v):
    return round(v / GRID) * GRID


def ssc_polygons(des: SSCDesign, x0=0.0, mirror=False):
    """Trapezoids (x along propagation).  mirror=True puts the facet at -x."""
    polys = []
    for (z0, z1, w0, w1) in des.teeth():
        a, b = (x0 + z0, x0 + z1) if not mirror else (x0 - z0, x0 - z1)
        pts = [(snap(a), snap(-w0 / 2)), (snap(b), snap(-w1 / 2)),
               (snap(b), snap(w1 / 2)), (snap(a), snap(w0 / 2))]
        polys.append(gdstk.Polygon(pts, *L_WG))
    return polys


def drc(polys):
    xs = sorted([(min(p.points[:, 0]), max(p.points[:, 0]),
                  min(np.ptp(p.points[p.points[:, 0] == p.points[:, 0].min(), 1]),
                      np.ptp(p.points[p.points[:, 0] == p.points[:, 0].max(), 1])))
                 for p in polys])
    tooth = min(b - a for a, b, _ in xs[1:])
    gaps = [xs[i + 1][0] - xs[i][1] for i in range(len(xs) - 1)]
    gaps = [g for g in gaps if g > 1e-6]
    width = min(w for _, _, w in xs)
    return dict(min_tooth_um=round(tooth, 4), min_gap_um=round(min(gaps), 4),
                min_width_um=round(width, 4), n_teeth=len(xs) - 1,
                length_um=round(xs[-1][1] - xs[0][0], 3))


def build(des: SSCDesign, out="ssc_2p5pct_3um_1550.gds"):
    lib = gdstk.Library("SSC_2P5PCT_3UM_1550", unit=1e-6, precision=1e-9)
    Ldev = des.teeth()[-1][1]

    c_ssc = lib.new_cell("SSC_2p5_3um_1550")
    polys = ssc_polygons(des)
    c_ssc.add(*polys)
    c_ssc.add(gdstk.rectangle((snap(Ldev), -15), (snap(Ldev) + 0.2, 15), *L_FACET))
    c_ssc.add(*gdstk.text("FACET", 4, (Ldev - 30, 12), layer=L_TEXT[0], datatype=L_TEXT[1]))
    c_ssc.add(*gdstk.text("CHIP WG 3um", 4, (0, 12), layer=L_TEXT[0], datatype=L_TEXT[1]))
    rep = drc(polys)

    Lmid = 3000.0
    c_b2b = lib.new_cell("SSC_B2B_TEST")
    left = Ldev                         # facet of mirrored SSC at x=0
    c_b2b.add(*ssc_polygons(des, x0=left, mirror=True))
    c_b2b.add(gdstk.rectangle((snap(left), -1.5), (snap(left + Lmid), 1.5), *L_WG))
    c_b2b.add(*ssc_polygons(des, x0=left + Lmid))
    tot = 2 * Ldev + Lmid
    for xf in (0.0, tot):
        c_b2b.add(gdstk.rectangle((snap(xf) - 0.1, -15), (snap(xf) + 0.1, 15), *L_FACET))
    c_b2b.add(*gdstk.text("SSC B2B  L=%.0fum" % tot, 8, (left + 100, 14), layer=L_TEXT[0], datatype=L_TEXT[1]))

    c_ref = lib.new_cell("REF_STRAIGHT")
    c_ref.add(gdstk.rectangle((0, -1.5), (snap(tot), 1.5), *L_WG))
    for xf in (0.0, tot):
        c_ref.add(gdstk.rectangle((snap(xf) - 0.1, -15), (snap(xf) + 0.1, 15), *L_FACET))
    c_ref.add(*gdstk.text("REF 3um WG  L=%.0fum" % tot, 8, (100, 14), layer=L_TEXT[0], datatype=L_TEXT[1]))

    # alternative: solid lateral inverse taper 3 -> 0.55 um over 2 mm (sqrt law), 20 um lead
    c_alt = lib.new_cell("SSC_ALT_INV_TAPER")
    n = 400
    s_ = (np.arange(n + 1) / n)
    xs = 20.0 + 2000.0 * s_
    ws = 3.0 + (0.55 - 3.0) * np.sqrt(s_)
    pts = [(0, -1.5)] + [(snap(x), snap(-w / 2)) for x, w in zip(xs, ws)] + \
          [(snap(x), snap(w / 2)) for x, w in zip(xs[::-1], ws[::-1])] + [(0, 1.5)]
    c_alt.add(gdstk.Polygon(pts, *L_WG))
    c_alt.add(gdstk.rectangle((snap(xs[-1]), -15), (snap(xs[-1]) + 0.2, 15), *L_FACET))
    c_alt.add(*gdstk.text("ALT solid inverse taper tip 0.55um", 6, (100, 12), layer=L_TEXT[0], datatype=L_TEXT[1]))

    top = lib.new_cell("TOP")
    top.add(gdstk.Reference(c_alt, (0, 750)))
    top.add(gdstk.Reference(c_ssc, (0, 500)))
    top.add(gdstk.Reference(c_b2b, (0, 250)))
    top.add(gdstk.Reference(c_ref, (0, 0)))
    lib.write_gds(os.path.join(HERE, out))
    return rep, tot


def plot_layout(des: SSCDesign, png="ssc_layout.png"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon as MP
    segs = des.teeth()
    Ldev = segs[-1][1]
    fig, axs = plt.subplots(3, 1, figsize=(11, 7.5), gridspec_kw=dict(height_ratios=[1.2, 1, 1]))
    windows = [(-10, Ldev + 10, "Full SSC (x stretched)"), (0, 60, "Chip side: solid 3 um lead -> first gaps (0.5 um)"),
               (Ldev - 40, Ldev + 4, "Facet side: last teeth (0.5 um) at SMF facet")]
    for ax, (a, b, ttl) in zip(axs, windows):
        for (z0, z1, w0, w1) in segs:
            if z1 < a or z0 > b:
                continue
            ax.add_patch(MP([(z0, -w0 / 2), (z1, -w1 / 2), (z1, w1 / 2), (z0, w0 / 2)],
                            closed=True, fc="#3b5bdb", ec="#1c2f80", lw=0.4))
        ax.axvline(Ldev, color="#e03131", lw=1, ls="--")
        ax.set_xlim(a, b); ax.set_ylim(-4, 4)
        ax.set_title(ttl, fontsize=10); ax.set_ylabel("y (um)")
    axs[-1].set_xlabel("x along propagation (um)   [red dashed: facet / SMF]")
    fig.tight_layout(); fig.savefig(os.path.join(HERE, png), dpi=140); plt.close(fig)


if __name__ == "__main__":
    des, src = load_design(sys.argv[1] if len(sys.argv) > 1 else None)
    rep, tot = build(des)
    plot_layout(des)
    print("design from", src, des.to_dict())
    print("DRC:", rep, " B2B length:", tot)
    json.dump(dict(design=des.to_dict(), drc=rep, b2b_length_um=tot),
              open(os.path.join(HERE, "gds_report.json"), "w"), indent=2)
