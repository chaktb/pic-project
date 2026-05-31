AWG WDM Design Utility V2.01 (Python/Tkinter Port)
====================================================
Original: (C) 2000 by Jang-Uk Shin and ETRI
Python port: 2024

REQUIREMENTS
------------
Python 3.8+
Tkinter (included in standard Python on macOS/Windows/Linux)
No additional packages needed.

HOW TO RUN
----------
    python3 main.py

or on Windows:
    python main.py

FILE STRUCTURE
--------------
main.py          - Application entry point
globals.py       - Global variables and constants
calc_slab.py     - Slab waveguide mode calculations (Calc_Slab.bas)
awg_calc.py      - AWG core calculations (Set_Globals.bas)

frm_variables.py - Main design parameter form (Frm_Variables.frm)
frm_index.py     - Refractive index calculator (frm_Index.frm)
frm_param.py     - Design parameters & simulation (Frm_Param.frm)
frm_plot.py      - Layout parameters (Frm_Plot.frm)
frm_layout.py    - Device layout drawing (Frm_Layout.frm)
frm_simulation.py - Spectral response plot (frm_Simulation.frm)
frm_cplcoeff.py  - Coupling coefficient plot (frm_CplCoeff.frm)
frm_add.py       - Taper/parabolic waveguide items (Frm_Add.frm)
frm_rindex.py    - Refractive index data (frm_RIndex.frm)
frm_coeff_edit.py - Per-waveguide parameter editor (frm_CoeffEdit.frm)
frm_about.py     - About dialog (Frm_ID.frm)

WORKFLOW
--------
1. Start the application: python3 main.py
2. Enter wavelength parameters (center wavelength, channel spacing, # channels)
3. Enter waveguide parameters (width, height, polarization)
4. Use "Waveguide" menu to open the refractive index calculator
   - Enter material names and Sellmeier coefficients
   - Click "Calculate" to get effective indices
   - Click "OK" to apply to main form
5. Enter AWG parameters (# array waveguides, spacings)
6. Click "Simulate" to:
   - Calculate FSR, path difference, focal length
   - Open the Design Parameters window
7. In Design Parameters window:
   - Choose coupling profile (Uniform, Gaussian, or BPM)
   - Set simulation wavelength range and points
   - Click "Simulate" to calculate spectral response
8. Click "Draw" to:
   - Open Layout Parameters window
   - Set device geometry parameters
   - Click "Draw" to render the device layout
9. From Layout window, use File > Make DXF to export

FILE FORMAT (.awg)
------------------
Plain text file, one parameter per line.
Compatible with original VB version .awg files.

PHYSICS
-------
- Effective Index Method (EIM) for waveguide mode calculation
- Hybrid EIM (HEIM) with mode table lookup available
- Bisection method for characteristic equation solving
- Transfer matrix approach for spectral simulation
- Rowland circle geometry for slab design
