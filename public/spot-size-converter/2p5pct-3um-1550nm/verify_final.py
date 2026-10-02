"""
verify_final.py -- fine-grid verification of the GA-optimized SSC
=================================================================
  1. 3D BPM through the REAL tooth/gap geometry at dx = 0.10 um, dz <= 0.2 um
     (converged grid), with x-z / y-z field records
  2. effective-medium BPM and local-mode EME (product of adjacent local-mode
     overlaps x facet overlap) for comparison
  3. facet-cut position tolerance (where the polish line lands inside the last period)
  4. SMF lateral / vertical misalignment tolerance
  5. wavelength sweep 1500-1600 nm
  6. lithography CD bias +-0.1 um (teeth longer/wider vs shorter/narrower)
Writes final_design.json, verify_results.json, field_records.npz
"""
import json, math, time, copy, os
import numpy as np
from ssc_core import Platform, Grid, BPM, SSCDesign, mode_solve, db, d4sigma, segmented_slices

HERE = os.path.dirname(os.path.abspath(__file__))
DX = float(os.environ.get("SSC_DXF", 0.10))


def load():
    for f in ("ga_stage3_refine.json", "ga_stage2_seg.json"):
        p = os.path.join(HERE, f)
        if os.path.exists(p):
            d = json.load(open(p))
            g = d["best_genes"]
            return SSCDesign(L_um=g["L_um"], pitch_um=g["pitch_um"], d_facet=g["d_facet"],
                             w_facet_um=g["w_facet_um"], p_d=g["p_d"], p_w=g["p_w"]), f, d
    raise FileNotFoundError


class Biased(SSCDesign):
    """CD bias b: every tooth grows by b in length (gaps shrink) and in width."""
    bias: float = 0.0

    def teeth(self, grid_snap=0.001):
        out = []
        for i, (z0, z1, w0, w1) in enumerate(super().teeth(grid_snap)):
            if i == 0:
                out.append((z0, z1 + self.bias / 2, w0 + self.bias, w1 + self.bias))
            else:
                out.append((z0 - self.bias / 2, z1 + self.bias / 2, w0 + self.bias, w1 + self.bias))
        return out


def run_case(des, plat, grid, record=0, cut=None):
    bpm = BPM(grid, plat, pml_um=3.0)
    smf = grid.gaussian(plat.smf_mfd_um / 2)
    ne, chip = mode_solve(grid, plat, grid.n2_rect(plat, plat.core_w_um))
    sl = segmented_slices(des, grid, plat, 0.2)
    if cut is not None:                   # extra cladding-free propagation or truncated last tooth
        sl = sl + cut
    psi, rec = bpm.run(chip.copy(), sl, record_every=record or None)
    return psi, rec, smf, chip


def main():
    t0 = time.time()
    plat = Platform()
    des, src, gares = load()
    grid = Grid(16.0, 14.0, DX)
    res = dict(source=src, design=des.to_dict(), drc=des.drc(), grid_dx_um=DX)

    # ---- 1. fine-grid true-segment BPM
    psi, rec, smf, chip = run_case(des, plat, grid, record=20)
    eta = grid.overlap_raw(psi, smf)
    res["seg_bpm_loss_dB"] = db(eta)
    res["seg_bpm_eta"] = eta
    res["power_at_facet"] = grid.power(psi)
    res["bare_loss_dB"] = db(grid.overlap(chip, smf))
    res["chip_d4s_um"] = d4sigma(grid, chip)
    zs, rx, ry = rec
    np.savez_compressed(os.path.join(HERE, "field_records.npz"), z=zs, x=grid.x, y=grid.y,
                        rx=rx.astype(np.float32), ry=ry.astype(np.float32),
                        facet=np.abs(psi).astype(np.float32))
    print(f"[1] seg BPM dx={DX}: {res['seg_bpm_loss_dB']:.4f} dB  (bare {res['bare_loss_dB']:.3f} dB)  "
          f"P_facet={res['power_at_facet']:.4f}  {time.time()-t0:.0f}s", flush=True)

    # ---- 2. effective medium BPM + local-mode EME
    bpm = BPM(grid, plat, pml_um=3.0)
    track = des.em_track(1.0)
    sl = []
    for zl, w, d in track:
        n = max(1, int(round(zl)))
        sl.append((grid.n2_rect(plat, w, d), zl / n, n))
    psi_em, _ = bpm.run(chip.copy(), sl)
    res["em_bpm_loss_dB"] = db(grid.overlap_raw(psi_em, smf))
    duty, width = des.periods()
    idx = np.unique(np.linspace(0, len(duty) - 1, 80).astype(int))
    modes = [chip] + [mode_solve(grid, plat, grid.n2_rect(plat, width[i], duty[i]))[1] for i in idx]
    ov = [grid.overlap(modes[i], modes[i + 1]) for i in range(len(modes) - 1)]
    nef, fac = mode_solve(grid, plat, grid.n2_rect(plat, width[-1], duty[-1]))
    res["eme"] = dict(junction_dB=db(ov[0]), taper_dB=db(np.prod(ov)),
                      facet_mode_smf_dB=db(grid.overlap(fac, smf)),
                      total_dB=db(np.prod(ov) * grid.overlap(fac, smf)),
                      facet_neff=nef, facet_margin=nef - plat.n_clad,
                      facet_d4s_um=d4sigma(grid, fac))
    print(f"[2] EM-BPM {res['em_bpm_loss_dB']:.4f} dB   EME {res['eme']['total_dB']:.4f} dB "
          f"(facet-mode {res['eme']['facet_mode_smf_dB']:.4f}, D4s {res['eme']['facet_d4s_um']})", flush=True)

    # ---- 3. facet cut position: the facet lands x um past the nominal cut
    #        x in (0, gap): polish line falls in the last gap -> add cladding propagation
    #        negative: cut back into the last tooth (shorter last tooth)
    cuts = {}
    n2c = np.full((grid.nx, grid.ny), plat.n_clad ** 2)
    P = des.pitch_um
    for off in (-0.25, 0.0, 0.5, 1.0, 1.5):
        if off < 0:                                   # last tooth trimmed by |off|
            segs = segmented_slices(des, grid, plat, 0.2)
            n2last, dzl, nl = segs[-1]
            keep = max(0, nl - int(round(-off / dzl)))
            segs = segs[:-1] + ([(n2last, dzl, keep)] if keep else [])
            psi_c, _ = BPM(grid, plat, pml_um=3.0).run(chip.copy(), segs)
        elif off == 0:
            psi_c = psi
        else:                                         # cut lands inside the following gap
            n = max(1, int(round(off / 0.1)))
            psi_c, _ = BPM(grid, plat, pml_um=3.0).run(psi.copy(), [(n2c, off / n, n)])
        cuts[f"{off:+.2f}"] = db(grid.overlap_raw(psi_c, smf))
    res["facet_cut_dB"] = cuts
    print(f"[3] facet cut offset (um) -> loss: {cuts}", flush=True)

    # ---- 4. SMF misalignment (shifted Gaussian)
    mis = {}
    for dxs in (0.5, 1.0, 1.5, 2.0):
        gsx = np.exp(-((grid.X - dxs) ** 2 + grid.Y ** 2) / (plat.smf_mfd_um / 2) ** 2)
        gsy = np.exp(-(grid.X ** 2 + (grid.Y - dxs) ** 2) / (plat.smf_mfd_um / 2) ** 2)
        mis[f"{dxs:.1f}"] = (db(grid.overlap_raw(psi, grid.normalize(gsx.astype(complex)))),
                             db(grid.overlap_raw(psi, grid.normalize(gsy.astype(complex)))))
    res["misalign_dB_xy"] = mis
    print(f"[4] misalignment (dx -> (loss_x, loss_y)): {mis}", flush=True)

    # ---- 5. wavelength sweep (SMF MFD scales ~ linearly near 1550: 10.4 um * (lam/1.55)^1.0 approx.)
    wl = {}
    for lam in (1.50, 1.525, 1.55, 1.575, 1.60):
        pl = Platform(lam_um=lam, smf_mfd_um=10.4 * (lam / 1.55) ** 1.0)
        if abs(lam - 1.55) < 1e-9:
            wl[f"{lam*1000:.0f}"] = res["seg_bpm_loss_dB"]; continue
        p2, _, s2, _ = run_case(des, pl, grid)
        wl[f"{lam*1000:.0f}"] = db(grid.overlap_raw(p2, s2))
        print(f"    lam {lam*1000:.0f} nm: {wl[f'{lam*1000:.0f}']:.4f} dB", flush=True)
    res["wavelength_dB"] = wl

    # ---- 6. CD bias (only shrinking bias keeps >= 0.5 um nominal rule -- reported both ways)
    cd = {}
    for b in (-0.1, -0.05, 0.05, 0.1):
        dd = Biased(**{k: v for k, v in dict(L_um=des.L_um, pitch_um=des.pitch_um, d_facet=des.d_facet,
                                              w_facet_um=des.w_facet_um, p_d=des.p_d, p_w=des.p_w,
                                              d_start=des.d_start).items()})
        dd.bias = b
        p3, _, s3, _ = run_case(dd, plat, grid)
        cd[f"{b:+.2f}"] = db(grid.overlap_raw(p3, s3))
        print(f"    CD bias {b:+.2f} um: {cd[f'{b:+.2f}']:.4f} dB", flush=True)
    res["cd_bias_dB"] = cd

    res["runtime_s"] = time.time() - t0
    json.dump(res, open(os.path.join(HERE, "verify_results.json"), "w"), indent=2, default=float)
    json.dump(dict(design=des.to_dict(), drc=des.drc(), loss_dB=res["seg_bpm_loss_dB"], source=src),
              open(os.path.join(HERE, "final_design.json"), "w"), indent=2, default=float)
    print("done", res["runtime_s"])


if __name__ == "__main__":
    main()
