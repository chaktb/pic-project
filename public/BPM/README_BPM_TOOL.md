# GDS BPM Tool (Tkinter) — 2D (EIM) and 3D (FD-ADI)

## What this tool does
- Opens a GDS file and extracts polygons from a selected layer/datatype.
- Assumes **GDS X = propagation axis z** and **GDS Y = transverse axis x**.
- Two simulation modes, selectable in the GUI:
  - **2D (EIM)** — the original mode. Builds a 2D effective-index map from the
    core polygons (slab TE0 effective index per lateral segment) and propagates a
    scalar field `ψ(x, z)` with a split-step Fourier BPM. The vertical (y)
    direction is collapsed by the effective-index method.
  - **3D (FD-ADI)** — a genuine three-dimensional solver. Propagates a scalar
    field `ψ(x, y, z)` over two transverse axes. The vertical (y) refractive
    index is defined by a **vertical layer stack** (see below). Uses a
    Peaceman–Rachford **ADI** scheme (unconditionally stable, 2nd order in z)
    with a **stretched-coordinate PML** absorbing boundary. Lives in `bpm3d.py`.
- Launches an input Gaussian mode (1D in 2D mode, 2D in 3D mode).
- Calculates total / per-port / target-channel output power and insertion loss.

## 3D mode and the vertical layer stack
A GDS file only contains the lateral (x–z) layout, so the out-of-plane (y)
structure must be supplied. In the **3D / Vertical Layer Stack** panel each entry
maps to a vertical slab `(y-bottom, thickness, index)`:
- **Patterned** layers take their lateral shape from a GDS layer (the same core
  mask the 2D tool builds). The index is applied only where that GDS layer's
  polygons exist, within the slab's y-range.
- **Non-patterned** layers (substrate, buffer, top cladding) fill the whole
  cross-section for their y-range.
Click **Auto** to generate a default substrate / core / top-clad stack from the
current indices and core thickness, then **Add / Edit / Remove** to refine it.
The paraxial reference index `n_ref` is estimated with a two-step effective
index (vertical slab, then lateral slab).

3D outputs: input index cross-section `n(x,y)`, top view `∫|ψ|²dy (z,x)`, side
view `∫|ψ|²dx (z,y)`, and the output cross-section `|ψ(x,y)|²`.

### Polarization — 3D (Scalar / TE / TM / Full-vector)
The 3D panel has a **Polarization** selector: `Scalar`, `TE`, `TM`, or `Full`.

- **TE / TM** use a **semi-vectorial** formulation — the transverse operator
  along the dominant E-field direction (TE → lateral x, TM → vertical y) carries
  the index-discontinuity correction `∂/∂u[(1/n²)∂(n²ψ)/∂u]`, so TE and TM modes
  get their distinct effective indices (the splitting grows with index contrast).
  In a uniform region the correction reduces to the scalar operator, so the PML
  is unaffected.
- **Full** is **full-vectorial with cross-polarization coupling**: it propagates
  both transverse components `(Ex, Ey)` together. The diagonal blocks are the
  semi-vectorial TE (for Ex) and TM (for Ey) operators; the off-diagonal coupling
  `Pxy = ∂/∂x[(∂ln n²/∂y)Ey]`, `Pyx = ∂/∂y[(∂ln n²/∂x)Ex]` is applied with
  Strang splitting. This coupling is non-zero only where both transverse index
  gradients meet — at **waveguide corners** — which is the physical origin of
  polarization conversion. Pick which component to launch with **Launch(full):
  Ex / Ey**; the other component starts empty and is filled by the coupling. The
  reported **polarization conversion** is the fraction of output power in the
  minor component, and the panel plots `|Ex(x,y)|²` and the converted `|Ey(x,y)|²`.
  (Axis-aligned waveguides barely convert; tilted/rotated cross-sections convert
  strongly.)

### Eigenmode input — 3D
The **Input field** selector chooses `Gaussian` or `Eigenmode`. With `Eigenmode`,
the tool first solves the **fundamental guided mode of the input cross-section**
by imaginary-distance BPM (reusing the ADI propagator — no extra dependency),
launches that mode, and reports its `n_eff`. A launched eigenmode keeps its shape
along a straight guide (minimal breathing), unlike a Gaussian. The mode solve is
polarization-aware (the reported `n_eff` differs for TE vs TM). For **Full**
polarization it solves the **full-vector hybrid mode** (coupled Ex/Ey via coupled
imaginary-distance iteration), so the launched field is a true vector eigenstate
of the chosen launch component — it propagates with stable polarization content
(no spurious conversion), whereas a non-eigen launch beats between polarizations.

### Output modal-loss matrix — 3D
Tick **Output eigenmode-overlap loss** to compute a **per-output-port modal-loss
matrix**. For each detected output port the tool isolates that port's guide (masks
the cross-section to a window around the port so the solver returns *that* port's
mode, not the global fundamental), solves its eigenmode, and couples the
propagated field into it — a physically meaningful **modal coupling loss** rather
than a window-integrated power. For **Full** polarization it solves both the TE
and TM mode of each port and reports two columns (Ex→TE, Ey→TM), so the matrix
shows how much co- and cross-polarized power each port receives. Costs one (or two,
for full-vector) extra mode solves per port.

### Wavelength sweep — 3D
With 3D selected, **Sweep Monitor Loss** runs a 3D BPM at each wavelength in the
Start/Stop/Step range and plots **total + selected-port insertion loss vs λ**, and
for Full polarization also **polarization-conversion (%) vs λ** on a second axis.
The sweep uses a Gaussian input. Tick **Material dispersion (silica)** to scale
every stack index per wavelength by the silica Sellmeier ratio
`n_silica(λ)/n_silica(λ_ref)` (a uniform approximation referenced to the main
Wavelength field); untick it to hold indices fixed (geometric dispersion only).
It is heavy — each wavelength is a full 3D run.

## Important limitations
This is a practical Python tool, not a full drop-in replacement for the original
ETRI C++ code. The 2D mode is **scalar**; the 3D mode offers scalar,
**semi-vectorial** TE/TM, and **full-vectorial** (Ex+Ey with cross-polarization
coupling). All modes are **paraxial**. The full-vectorial cross term is treated
explicitly (Strang-split), which is accurate for the usual corner-driven coupling
but not a fully-implicit vector scheme; eigenmode input for full-vector launches
the semi-vectorial mode of the chosen component (a full-vector hybrid-mode solver
is a possible future addition). The 3D ADI solver uses a Python-loop Thomas solve
per z-step, so very large grids are slower than a compiled FD-BPM (TE/TM/full
rebuild the index-dependent operators each z-step).

Good for: directional trend checking, channel routing / loss comparison, and
layout-to-BPM experiments with a real vertical stack.

## Install
```bash
pip install numpy matplotlib shapely gdstk
```
(No extra dependency for 3D — `bpm3d.py` uses only numpy; the Thomas solver is
pure numpy, no scipy.)

## Run
```bash
python bpm_gui_rev08.py
```

## Test (headless, numpy only)
```bash
python test_bpm3d.py   # free-space diffraction, straight wg, PML, GUI data pipeline
```

## Recommended GDS usage
- Put the **waveguide core outline** on one dedicated layer.
- Keep propagation mainly along the GDS **X direction**.
- Use the GDS **Y direction** as the transverse direction.
- For multi-channel output, set `Output center` to the channel center you want to analyze.

## Suggested next upgrades
1. Auto-detect input/output port centers from the GDS edges.
2. Replace Gaussian trial mode with an eigenmode solver.
3. Add TE/TM semi-vectorial BPM.
4. Add bend loss correction and effective-index curvature model.
5. Support multi-layer stack and different material regions.
6. Export field snapshots and channel tables automatically.
