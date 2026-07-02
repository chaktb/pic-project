import numpy as np, pickle
from ssc_ga_optimizer import Platform, GridConfig
from seg_ssc import SegGABounds, SegSSCGeneticOptimizer, evaluate_seg_gene, seg_gene_to_gds

plat = Platform(lam_um=1.301, delta=0.015, core_thick_um=5.0, smf_mfd_um=9.2)
grid = GridConfig(dx_um=0.25, dy_um=0.25, dz_um=0.6, x_half_um=16, y_half_um=13, pml_um=3.0, save_every=20)
bounds = SegGABounds(pitch_um=(3.0,8.0), duty_start=(0.75,0.98), duty_end=(0.10,0.40),
                     w_start=(5.0,5.0), w_end=8.5, n_seg=(35,100),
                     profiles=("linear","quad","cos"))
print('1.5%% 5x5um: n_clad=%.5f n_core=%.5f'%(plat.n_clad,plat.n_core))
ga = SegSSCGeneticOptimizer(plat, grid, bounds, pop_size=12, n_gen=8, seed=7, leadins=[(50.0,5.0)])
r = ga.run()
g = r['best_gene']
print('=== BEST ===')
print('loss=%.4f dB eff=%.3f%%'%(r['best_loss_dB'], r['best_eval']['efficiency']*100))
print('pitch=%.3f nseg=%d duty %.3f->%.3f(%s) width %.1f->%.1f(%s) L=%.1f'%(
    g.pitch_um,g.n_seg,g.duty_start,g.duty_end,g.duty_profile,g.w_start,g.w_end,g.width_profile,g.total_length))
pickle.dump((g,r),open('best15.pkl','wb'))
seg_gene_to_gds(g,'ssc_optimized.gds',y_offset_um=200.0)
print('saved gds')
