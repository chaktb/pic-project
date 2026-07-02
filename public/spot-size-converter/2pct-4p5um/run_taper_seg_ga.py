"""GA for tapered-segment SSC: start 4.5um solid chip -> end 8.5um SMF facet.
Fixed: w_start=4.5, w_end=8.5, thickness 4.5um, single 4.5um solid lead-in.
GA free params: pitch, n_seg, duty_start, duty_end, duty/width ramp profile.
"""
import numpy as np
from ssc_ga_optimizer import Platform, GridConfig
from seg_ssc import (SegSSCGene, SegGABounds, SegSSCGeneticOptimizer,
                     evaluate_seg_gene, seg_gene_to_gds)

def main(pop, gen, seed, fast, gds_out):
    plat = Platform()
    if fast:
        grid = GridConfig(dx_um=0.25, dy_um=0.25, dz_um=0.6,
                          x_half_um=16, y_half_um=13, pml_um=3.0, save_every=20)
    else:
        grid = GridConfig(dx_um=0.18, dy_um=0.18, dz_um=0.4,
                          x_half_um=18, y_half_um=15, pml_um=3.5, save_every=15)

    # bounds: start width fixed at 4.5 (chip), bus fixed 8.5 (SMF facet)
    bounds = SegGABounds(
        pitch_um=(3.0, 8.0),
        duty_start=(0.75, 0.98),
        duty_end=(0.10, 0.40),
        w_start=(4.5, 4.5),        # FIXED chip width
        w_end=8.5,                 # FIXED SMF-facet width
        n_seg=(35, 100),
        profiles=("linear", "quad", "cos"),
    )
    print(f"Platform: n_clad={plat.n_clad:.5f} n_core={plat.n_core:.5f} "
          f"t={plat.core_thick_um}um SMF w0={plat.smf_w0_um:.2f}um")
    print("Fixed: w_start=4.5um (chip, solid)  w_end=8.5um (SMF facet)  thickness=4.5um\n")

    ga = SegSSCGeneticOptimizer(
        plat, grid, bounds, pop_size=pop, n_gen=gen, seed=seed,
        leadins=[(50.0, 4.5)],     # single solid 4.5um chip lead-in
    )
    res = ga.run()
    g = res["best_gene"]
    print("\n=== BEST TAPERED-SEGMENT SSC ===")
    print(f"coupling loss @1301nm : {res['best_loss_dB']:.4f} dB")
    print(f"efficiency            : {res['best_eval']['efficiency']*100:.3f} %")
    print(f"pitch / n_seg         : {g.pitch_um:.3f} um / {g.n_seg}")
    print(f"duty  start->end      : {g.duty_start:.3f} -> {g.duty_end:.3f} ({g.duty_profile})")
    print(f"width start->end      : {g.w_start:.3f} -> {g.w_end:.3f} um ({g.width_profile})")
    print(f"total length          : {g.total_length:.1f} um")
    seg_gene_to_gds(g, gds_out, y_offset_um=200.0)
    print(f"GDS written           : {gds_out}")
    return plat, grid, g, res

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--pop", type=int, default=14)
    ap.add_argument("--gen", type=int, default=10)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--fast", action="store_true")
    ap.add_argument("--gds", default="taper_seg_best.gds")
    a = ap.parse_args()
    main(a.pop, a.gen, a.seed, a.fast, a.gds)
