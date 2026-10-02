"""
verify_widen.py -- fine-grid (dx = 0.10 um) verification of the GA-optimized
widening, minimum-pitch segmented SSC (seg_vp.WideningSSC).

  1. true tooth/gap 3D BPM, x-z / y-z field records, facet field
  2. facet-mode / SMF overlap (duty-averaged cross-section at the facet)
  3. facet cut position, SMF misalignment, wavelength 1500-1600 nm, CD bias +-0.1 um
  4. reference: previous uniform-width (3 um) duty-only design, same grid
Writes final_design.json, verify_results.json, field_records.npz
"""
import json, os, time
import numpy as np
from ssc_core import Platform, Grid, BPM, SSCDesign, mode_solve, db, d4sigma, segmented_slices
from seg_vp import WideningSSC, from_genes, NAMES

H = os.path.dirname(os.path.abspath(__file__))
DX = float(os.environ.get("SSC_DXF", 0.10))


class Biased(WideningSSC):
    bias: float = 0.0

    def teeth(self, grid_snap=0.001):
        out = []
        for i, (z0, z1, w0, w1) in enumerate(super().teeth(grid_snap)):
            b = self.bias
            out.append((z0 if i == 0 else z0 - b / 2, z1 + b / 2, w0 + b, w1 + b))
        return out


def run(des, plat, grid, record=0, extra=None):
    bpm = BPM(grid, plat, pml_um=3.0)
    smf = grid.gaussian(plat.smf_mfd_um / 2)
    _, chip = mode_solve(grid, plat, grid.n2_rect(plat, plat.core_w_um))
    sl = segmented_slices(des, grid, plat, 0.2)
    psi, rec = bpm.run(chip.copy(), sl, record_every=record or None)
    return psi, rec, smf, chip


def main():
    t0 = time.time()
    plat = Platform(); grid = Grid(16.0, 14.0, DX)
    src = "ga_widen_B.json" if os.path.exists(os.path.join(H, "ga_widen_B.json")) else "ga_widen_A.json"
    G = json.load(open(os.path.join(H, src)))
    v = [G["best_genes"][k] for k in NAMES]
    des = from_genes(v)
    R = dict(source=src, genes=G["best_genes"], design=des.to_dict(), drc=des.drc(), grid_dx_um=DX)

    psi, rec, smf, chip = run(des, plat, grid, record=20)
    R["loss_dB"] = db(grid.overlap_raw(psi, smf)); R["eta"] = grid.overlap_raw(psi, smf)
    R["power_at_facet"] = grid.power(psi)
    R["bare_loss_dB"] = db(grid.overlap(chip, smf)); R["chip_d4s_um"] = d4sigma(grid, chip)
    R["out_d4s_um"] = d4sigma(grid, psi)
    zs, rx, ry = rec
    np.savez_compressed(os.path.join(H, "field_records.npz"), z=zs, x=grid.x, y=grid.y,
                        rx=rx.astype(np.float32), ry=ry.astype(np.float32), facet=np.abs(psi).astype(np.float32))
    print(f"[1] loss {R['loss_dB']:.4f} dB (bare {R['bare_loss_dB']:.3f})  P_facet {R['power_at_facet']:.4f}  "
          f"out D4s {R['out_d4s_um']}  {time.time()-t0:.0f}s", flush=True)

    nef, fac = mode_solve(grid, plat, grid.n2_rect(plat, des.w_facet_um, des.d_facet))
    R["facet_mode"] = dict(neff=nef, margin=nef - plat.n_clad, d4s_um=d4sigma(grid, fac),
                           smf_loss_dB=db(grid.overlap(fac, smf)))
    print(f"[2] facet mode: {R['facet_mode']}", flush=True)

    n2c = np.full((grid.nx, grid.ny), plat.n_clad ** 2); cuts = {"+0.00": R["loss_dB"]}
    for off in (0.25, 0.5, 1.0, 2.0):
        n = max(1, int(round(off / 0.1)))
        p2, _ = BPM(grid, plat, pml_um=3.0).run(psi.copy(), [(n2c, off / n, n)])
        cuts[f"+{off:.2f}"] = db(grid.overlap_raw(p2, smf))
    segs = segmented_slices(des, grid, plat, 0.2)
    n2l, dzl, nl = segs[-1]
    keep = max(0, nl - int(round(0.25 / dzl)))
    p2, _ = BPM(grid, plat, pml_um=3.0).run(chip.copy(), segs[:-1] + ([(n2l, dzl, keep)] if keep else []))
    cuts["-0.25"] = db(grid.overlap_raw(p2, smf))
    R["facet_cut_dB"] = cuts; print(f"[3] facet cut: {cuts}", flush=True)

    mis = {}
    for d in (0.5, 1.0, 1.5, 2.0):
        gx = grid.normalize(np.exp(-((grid.X - d) ** 2 + grid.Y ** 2) / 5.2 ** 2).astype(complex))
        gy = grid.normalize(np.exp(-(grid.X ** 2 + (grid.Y - d) ** 2) / 5.2 ** 2).astype(complex))
        mis[f"{d:.1f}"] = (db(grid.overlap_raw(psi, gx)), db(grid.overlap_raw(psi, gy)))
    R["misalign_dB_xy"] = mis; print(f"[4] misalign: {mis}", flush=True)

    wl = {"1550": R["loss_dB"]}
    for lam in (1.50, 1.525, 1.575, 1.60):
        pl = Platform(lam_um=lam, smf_mfd_um=10.4 * lam / 1.55)
        p3, _, s3, _ = run(des, pl, grid)
        wl[f"{lam*1000:.0f}"] = db(grid.overlap_raw(p3, s3)); print(f"    {lam*1000:.0f} nm {wl[f'{lam*1000:.0f}']:.4f}", flush=True)
    R["wavelength_dB"] = dict(sorted(wl.items()))

    cd = {}
    for b in (-0.1, -0.05, 0.05, 0.1):
        dd = Biased(**des.to_dict_kwargs()) if hasattr(des, "to_dict_kwargs") else Biased(
            L_um=des.L_um, d_start=des.d_start, d_facet=des.d_facet, w_facet_um=des.w_facet_um,
            p_d=des.p_d, p_w=des.p_w, P_min_um=des.P_min_um)
        dd.bias = b
        p4, _, s4, _ = run(dd, plat, grid)
        cd[f"{b:+.2f}"] = db(grid.overlap_raw(p4, s4)); print(f"    CD {b:+.2f}: {cd[f'{b:+.2f}']:.4f}", flush=True)
    R["cd_bias_dB"] = cd

    # reference: previous uniform-width duty-only design
    old = json.load(open(os.path.join(H, "archive_uniform", "final_design_uniform.json")))["design"]
    od = SSCDesign(L_um=old["n_periods"] * old["pitch_um"], pitch_um=old["pitch_um"], d_facet=old["d_facet"],
                   w_facet_um=old["w_facet_um"], p_d=old["p_d"], p_w=old["p_w"], d_start=old["d_start"])
    p5, _, s5, _ = run(od, plat, grid)
    R["reference_uniform_dB"] = db(grid.overlap_raw(p5, s5))
    print(f"[ref] uniform-width design: {R['reference_uniform_dB']:.4f} dB", flush=True)

    R["runtime_s"] = time.time() - t0
    json.dump(R, open(os.path.join(H, "verify_results.json"), "w"), indent=2, default=float)
    json.dump(dict(design=des.to_dict(), genes=G["best_genes"], drc=des.drc(), loss_dB=R["loss_dB"]),
              open(os.path.join(H, "final_design.json"), "w"), indent=2, default=float)
    print("done", R["runtime_s"])


if __name__ == "__main__":
    main()
