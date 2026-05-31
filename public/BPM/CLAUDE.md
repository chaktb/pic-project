# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Scope

This directory is the **GDS BPM Tool** — a Tkinter desktop app that loads a GDS layout and runs a scalar beam-propagation simulation over it, in either a **2D (effective-index)** or **3D** mode. It lives inside the larger `pic-project` photonic-circuits repo but is self-contained; everything it needs is here.

- `bpm_gui_rev08.py` — the GUI application + the original 2D solver (`BPM1DSolver`).
- `bpm3d.py` — companion module: the 3D solver (`BPM3DSolver`, FD-ADI + PML) and the vertical layer-stack model (`StackLayer`, `VerticalStack`, `default_silica_stack`). Pure numpy, no scipy.
- `test_bpm3d.py` — headless numerical tests (numpy only, no GUI/GDS): free-space diffraction vs analytic, straight-waveguide guiding + power conservation, PML absorption, semi-vectorial TE/TM `n_eff` vs analytic slab, eigenmode low-breathing, full-vector cross-pol coupling (separable→zero, power conservation, rotated-core conversion), and the GUI's 3D data pipeline.
- `mmi_2x2.gds` — sample input (a 2×2 MMI coupler layout) for manual testing.
- `README_BPM_TOOL.md` — user-facing description, limitations, and the GDS authoring conventions; read it for the physics caveats.

## Run

Use the project's Nix Python interpreter, and install the deps the solver/loader need (numpy + matplotlib ship with the Nix env; `shapely` and `gdstk` do not):

```bash
pip install shapely gdstk
~/.nix-profile/bin/python bpm_gui_rev08.py
```

`gdstk` is loaded lazily — `GDSGeometryLoader.available` reports whether it imported, and the UI shows `gdstk: OK/missing`. The app still launches without it; GDS loading just fails with an install hint. There is no CLI mode and no build step — GUI verification is launching it and loading a `.gds`.

The **3D solver and its data pipeline are unit-tested headlessly** (numpy only, no GUI/GDS/shapely):

```bash
python test_bpm3d.py
```

Prefer this for verifying solver changes — it's fast and needs no display. Note: in this sandbox the Nix env at `~/.nix-profile/bin/python` may be missing; a numpy-only interpreter (e.g. one of `/nix/store/*-python3-*-env/bin/python3`) is enough to run `test_bpm3d.py` since `bpm3d.py` imports only numpy.

## Architecture

The file is organized top-to-bottom as **data containers → GDS loader → BPM solver → Tkinter GUI**.

**GDS geometry (`GDSGeometryLoader`, `GDSLayerGeometry`, `GDSScanInfo`, `PortCandidate`).** `scan()` enumerates `(layer, datatype)` specs and their polygon counts so the user can pick which layer is the waveguide core; `load()` extracts that layer's polygons. `GDSLayerGeometry.__post_init__` **pre-buffers Shapely polygons once** (a perf fix — see the rev header) so per-z-slice point-in-polygon tests don't re-`.buffer()`.

**Axis convention (critical).** GDS **X = propagation axis z**, GDS **Y = transverse axis x**. The simulation marches along GDS-X. Keep this straight when reading any `z_val` / transverse-`x` code — they do not match the GDS coordinate names.

**Solver (`BPM1DSolver`).** Builds a 2D effective-index map by slicing the geometry along z, finding core segments at each slice (`_exact_segments_at_z` / `sample_core_mask`), and assigning each segment a slab TE0 effective index. The slab index comes from a bisection root-finder on the waveguide characteristic equation, memoized via the module-level `_cached_slab_te0_neff` (`lru_cache` keyed on rounded inputs; widths are rounded to ~4 nm so identical waveguide segments share a cache entry). Propagation is a scalar split-step BPM in `run()`, which accumulates `total_power` only every `power_stride` steps to bound memory. Index is derived from `n_clad` + `delta_percent` (`n_core_from_delta_percent`) or silica Sellmeier (`silica_sellmeier_n`).

**3D solver (`bpm3d.py`).** `BPM3DSolver` propagates `ψ(x,y,z)` with the paraxial equation using the **same field convention** as `BPM1DSolver` (`E = ψ·exp(+iβ_ref·z)`). The core is `_adi_step` — one Peaceman–Rachford **ADI** step (implicit x / explicit y, then the reverse, index potential `V` split symmetrically), shared by both forward propagation (`run`) and mode solving. Transverse operators come from `_transverse_ops`: **scalar** uses plain stretched-coordinate-PML Laplacians (`_laplacian_tridiag`, built once); **TE/TM** are *semi-vectorial* — the operator along the dominant-E direction (TE→x, TM→y) carries the index-discontinuity correction `∂/∂u[(1/n²)∂(n²ψ)/∂u]` (`_corrected_laplacian_tridiag`, rebuilt per z-slice because it depends on `n`). The correction reduces to the scalar operator in uniform regions, so it composes cleanly with the PML. Tridiag coefficients may be 1D (line-constant) or 2D (index-dependent); `thomas_batch` / `_thomas_axis` / `_tridiag_matvec` handle both. To bound memory the 3D index `n(x,y)` is built **per z-slice on the fly** via an `index_fn(iz)` closure — never a full `n(x,y,z)` array.

**Eigenmode solver (`solve_mode`).** Finds the fundamental guided mode of a fixed cross-section by **imaginary-distance ADI** — the same `_adi_step` with `dz = -1j·dξ`, iterated and renormalised so the highest-`n_eff` mode dominates; `n_eff` comes from a Rayleigh quotient using real, PML-free operators (`_transverse_ops_real`). It forces PML off internally (the complex coordinate stretch destabilises the imaginary-distance iteration; a guided mode is localised so PML is unneeded). Polarization-aware (TE/TM/scalar). Used to launch a real waveguide mode instead of a Gaussian.

**Full-vectorial (`run_fullvec`, `pol='full'`).** Propagates `(Ex, Ey)` together. The diagonal blocks are the semi-vectorial TE operator for Ex and TM for Ey (`_fullvec_diag_ops`); the off-diagonal cross-coupling `Pxy = ∂/∂x[(∂ln n²/∂y)Ey]`, `Pyx = ∂/∂y[(∂ln n²/∂x)Ex]` (`_cross_sources`, central differences) is applied explicitly with **Strang splitting** (cross ½-step · diagonal ADI · cross ½-step, factored into `_fullvec_step`). Coupling is non-zero only at corners — physical polarization conversion. `BPM3DResult` gains `psi_out_minor` (Ey) and `conversion` (minor/total power). **`solve_mode_fullvec`** finds the true hybrid (Ex,Ey) eigenmode by *coupled* imaginary-distance iteration (`_fullvec_step` with `dz=-1j·dξ`, combined-norm renormalisation, neff via `_rayleigh_neff_fullvec` which includes the cross terms).

**Output modal-loss matrix.** After propagation the GUI optionally builds a per-output-port matrix: for each port it *isolates* the guide (masks the output cross-section index to a lateral window around the port, so `solve_mode` returns that port's mode rather than the global fundamental), solves it, and couples the field in via `overlap_power`. Full-vector solves each port's TE and TM mode (temp `pol='TE'`/`'TM'` solvers) → columns (Ex→TE, Ey→TM), exposing cross-pol routing. Stored on `modal_matrix3d` and in JSON.

**3D wavelength sweep dispersion.** `run_wavelength_sweep_3d` optionally applies silica material dispersion: per wavelength it scales every stack index by `n_silica(λ)/n_silica(λ_ref)` (`scale_stack`, λ_ref = the main Wavelength field) and rebuilds `bg_profile`/`bands`/`n_ref`; masks are wavelength-independent so they're built once. n_ref reference width is snapshotted on the main thread (`_estimate_n_ref_3d(..., ref_w_um=...)`) to avoid touching tk vars from the worker.

**Vertical stack.** A `VerticalStack` of `StackLayer`s: *patterned* layers take their lateral shape from a GDS layer's core mask within a y-band; *non-patterned* layers fill the cross-section (substrate/cladding). `background_profile`, `patterned_bands`, and `make_y_grid` feed the `index_fn`.

**GUI (`BPMGui(tk.Tk)`).** `_build_ui` composes the panels via `_section_file / _section_bpm / _section_stack / _section_sweep / _section_ports / _section_actions`. The pipeline a user drives is: **open GDS → detect/select ports → (pick 2D or 3D mode; for 3D edit the vertical stack, polarization, input field) → run BPM → plot/save JSON**. `run_bpm` dispatches to `run_bpm_3d` when `sim_mode == "3D"` (which itself branches to `run_fullvec` when `pol == "full"`); `run_wavelength_sweep` dispatches to `run_wavelength_sweep_3d` (per-wavelength 3D run → IL + conversion vs λ, `_show_sweep3d_plot`). The 3D path reuses `BPM1DSolver.build_core_mask_grid` (which now accepts a precomputed `x_in`/`z_in`) to build a per-GDS-layer lateral mask on one shared grid, then assembles the `index_fn` and runs `BPM3DSolver`. Long operations follow a fixed `_run_async(work_fn, done_fn, ...)` pattern: `work` runs on a worker thread and reports progress through a callback; results come back through a queue drained by `_poll_worker_queue` on the Tk main loop, then `done` updates widgets. **Never touch Tk widgets from inside a `work` function** — only from `done` / the main thread (the 3D `work` reads all tk vars up front for this reason).

**Conventions.** The module docstring is a running changelog of bug-fix + performance revisions (rev04 header, file named rev08); when you fix something non-trivial, append a numbered entry there in the same style. Internal lengths are in micrometres (`_um` suffix throughout). Some comments are in Korean — that's expected, leave them.

## Note on the parent repo's auto-commit

The repo root has a Stop hook (`../../.claude/auto-build.py`) that runs after every Claude turn and **auto-commits and pushes the whole tree to `main`**. Anything you change or leave in this directory will be committed without a separate prompt. This directory has no `index.html`, so the hook's landing-page card generator ignores it.
