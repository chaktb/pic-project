"""
AWG Slab Waveguide Mode Calculations
Ported from Calc_Slab.bas (VB 5.0/6.0) to Python
"""

import math
from globals import PI, Max, Min, arcsin, nr_table, nr_index

# Physical constants
eps0 = 8.854187818e-12   # permittivity of free space [F/m]
mu0  = 1.25663706e-6     # permeability of free space [H/m]


# ─────────────────────────────────────────────────────────────
# Bisection root-finder for the waveguide characteristic equation
# ─────────────────────────────────────────────────────────────

def _charac(b, v, mode, n1, n2, n3, m):
    """
    Waveguide characteristic equation.
    b   : normalised propagation constant  (0..1)
    v   : V-number
    mode: 'TE' or 'TM'
    n1  : lower cladding index
    n2  : core index
    n3  : upper cladding index
    m   : mode order (0 for fundamental)
    """
    a = (n3 * n3 - n1 * n1) / (n2 * n2 - n3 * n3)

    if mode == "TE":
        f_val = (v * math.sqrt(1.0 - b)
                 - m * PI
                 - math.atan(math.sqrt((b + a) / (1.0 - b)))
                 - math.atan(math.sqrt(b / (1.0 - b))))
    elif mode == "TM":
        f_val = (v * math.sqrt(1.0 - b)
                 - m * PI
                 - math.atan((n2 * n2) / (n1 * n1) * math.sqrt((b + a) / (1.0 - b)))
                 - math.atan((n2 * n2) / (n3 * n3) * math.sqrt(b / (1.0 - b))))
    else:
        f_val = 0.0

    return f_val


def _bisec(v, mode, n1, n2, n3, order):
    """
    Find normalised propagation constant b using bisection method.
    Returns b in [0, 1).
    """
    ESP = 1e-10

    xL = 0.0
    tL = _charac(xL, v, mode, n1, n2, n3, order)

    xR = 1.0 - ESP
    tR = _charac(xR, v, mode, n1, n2, n3, order)

    if tL * tR > 0.0:
        return xL  # no root found – below cutoff

    while abs(xR - xL) > ESP:
        x = (xL + xR) / 2.0
        t = _charac(x, v, mode, n1, n2, n3, order)
        if t != 0.0:
            if tL * t < 0.0:
                xR = x
                tR = t
            elif tR * t < 0.0:
                xL = x
                tL = t

    return (xL + xR) / 2.0


# ─────────────────────────────────────────────────────────────
# Normal EIM (Effective Index Method)
# ─────────────────────────────────────────────────────────────

def Calc_NEIM(wl, pol, n_core, n_over, n_buff, wgd_h, wgd_w, buf_h,
              n_eff_out, n_side_out, betta_out):
    """
    Normal Effective Index Method.
    Returns (n_eff, n_side, betta).
    wl, wgd_h, wgd_w in SI units [m].
    """
    ko = 2.0 * PI / wl

    if pol == "X":
        mode, othermode = "TE", "TM"
    elif pol == "Y":
        mode, othermode = "TM", "TE"
    elif pol == "S":
        mode, othermode = "TE", "TE"
    else:
        raise ValueError(f"Unexpected polarization: {pol}")

    n_side = n_over

    # Vertical slab
    v = ko * wgd_h * math.sqrt(n_core ** 2 - Max(n_over, n_buff) ** 2)
    b = _bisec(v, mode, Min(n_over, n_buff), n_core, Max(n_over, n_buff), 0)
    n_eff = math.sqrt(b * (n_core ** 2 - Max(n_over, n_buff) ** 2)
                      + Max(n_over, n_buff) ** 2)

    # Horizontal slab
    v = ko * wgd_w * math.sqrt(n_eff ** 2 - n_side ** 2)
    b = _bisec(v, othermode, n_side, n_eff, n_side, 0)
    n_mode = math.sqrt(b * (n_eff ** 2 - n_side ** 2) + n_side ** 2)

    betta = n_mode * ko
    return n_mode, n_side, betta


# ─────────────────────────────────────────────────────────────
# Hybrid EIM (iterative)
# ─────────────────────────────────────────────────────────────

def Calc_HEIM(wl, pol, n_core, n_over, n_buff, wgd_h, wgd_w, buf_h,
              n_eff_out, n_side_out, betta_out, n_mode_init):
    """
    Hybrid Effective Index Method (iterative convergence).
    n_mode_init: initial guess for mode index (from table lookup).
    Returns (n_eff, n_side, betta).
    """
    n_err = 1e-10
    ko = 2.0 * PI / wl

    if pol == "X":
        mode, othermode = "TE", "TM"
    elif pol == "Y":
        mode, othermode = "TM", "TE"
    elif pol == "S":
        mode, othermode = "TE", "TE"
    else:
        raise ValueError(f"Unexpected polarization: {pol}")

    n_side = n_over

    # Vertical slab
    v = ko * wgd_h * math.sqrt(n_core ** 2 - Max(n_over, n_buff) ** 2)
    b = _bisec(v, mode, Min(n_over, n_buff), n_core, Max(n_over, n_buff), 0)
    n_eff = math.sqrt(b * (n_core ** 2 - Max(n_over, n_buff) ** 2)
                      + Max(n_over, n_buff) ** 2)

    # Horizontal slab (initial EIM)
    v = ko * wgd_w * math.sqrt(n_eff ** 2 - n_side ** 2)
    b = _bisec(v, othermode, n_side, n_eff, n_side, 0)
    n_eim = math.sqrt(b * (n_eff ** 2 - n_side ** 2) + n_side ** 2)

    n_mode = n_mode_init

    # Iterative adjustment
    max_iter = 1000
    for _ in range(max_iter):
        if abs(n_mode - n_eim) <= n_err:
            break
        n_side = n_side + 10.0 * (n_mode - n_eim)
        if n_side <= 0 or n_side >= n_eff:
            break
        v = ko * wgd_w * math.sqrt(n_eff ** 2 - n_side ** 2)
        b = _bisec(v, othermode, n_side, n_eff, n_side, 0)
        n_eim = math.sqrt(b * (n_eff ** 2 - n_side ** 2) + n_side ** 2)

    betta = n_eim * ko
    return n_eim, n_side, betta


# ─────────────────────────────────────────────────────────────
# Refractive index table interpolation (quadratic)
# ─────────────────────────────────────────────────────────────

def Calc_IndxTable(i_layer, wl_nm):
    """
    Quadratic interpolation of refractive index from lookup table.
    i_layer : layer index (0-based here; VB used 1-based)
    wl_nm   : wavelength in nm
    Returns  : interpolated refractive index
    """
    import globals as g

    # Convert VB 1-based layer to 0-based
    layer = i_layer   # already 0-based from Python callers

    n = g.nr_index[layer]
    if n <= 0:
        return 0.0

    tbl_wl = g.nr_table[layer][0]   # wavelength list (nm)
    tbl_n  = g.nr_table[layer][1]   # index list

    if wl_nm <= tbl_wl[0]:
        return tbl_n[0]
    if wl_nm >= tbl_wl[n - 1]:
        return tbl_n[n - 1]
    if n == 1:
        return tbl_n[0]

    # Find bracket
    i = 0
    while i < n and wl_nm > tbl_wl[i]:
        i += 1
    i -= 1  # now tbl_wl[i] <= wl_nm < tbl_wl[i+1]

    # Choose 3 points for quadratic interpolation
    if i == 0:
        i1, i2, i3 = 0, 1, 2
    elif i == n - 2:
        i1, i2, i3 = n - 3, n - 2, n - 1
    else:
        if (wl_nm - tbl_wl[i - 1]) < (tbl_wl[i + 1] - wl_nm):
            i1, i2, i3 = i - 1, i, i + 1
        else:
            i1, i2, i3 = i, i + 1, i + 2

    # Clamp
    if i3 >= n:
        i3 = n - 1
    if i1 < 0:
        i1 = 0

    x  = wl_nm
    x1 = tbl_wl[i1]
    x2 = tbl_wl[i2]
    x3 = tbl_wl[i3]
    y1 = tbl_n[i1]
    y2 = tbl_n[i2]
    y3 = tbl_n[i3]

    if x3 == x1:
        return y2

    result = (y1 * (x - x3) * (2 * x - x1 - x3)
              + 4 * y2 * (x - x1) * (x3 - x)
              + y3 * (x - x1) * (2 * x - x1 - x3)) / ((x3 - x1) ** 2)
    return result


# ─────────────────────────────────────────────────────────────
# Mode profile calculation and save
# ─────────────────────────────────────────────────────────────

def Save_ModePrf(filepath, wl, pol, d, n_core_eff, n_side):
    """
    Calculate 1-D mode profile and save to file.
    wl, d in SI units [m].
    """
    c_light = 1.0 / math.sqrt(eps0 * mu0)

    if pol == "X":
        mode = "TM"
    elif pol == "Y":
        mode = "TE"
    elif pol == "S":
        mode = "TE"
    else:
        raise ValueError(f"Unexpected polarization: {pol}")

    n2 = n_core_eff
    n3 = n_side
    n1 = n_side

    ko    = 2.0 * PI / wl
    omega = ko * c_light

    v = ko * d * math.sqrt(n2 ** 2 - n3 ** 2)
    b = _bisec(v, mode, n1, n2, n3, 0)

    if b == 0:
        raise ValueError("Mode profile not exist – mode cut off")

    n_beta = math.sqrt(b * (n2 ** 2 - n3 ** 2) + n3 ** 2)
    p    = ko * math.sqrt(n_beta ** 2 - n1 ** 2)
    q    = ko * math.sqrt(n_beta ** 2 - n3 ** 2)
    h    = ko * math.sqrt(n2 ** 2 - n_beta ** 2)
    beta = ko * n_beta

    dx = d / 25.0
    XX  = []
    Prf = []

    if mode == "TE":
        c_amp = 2 * h * math.sqrt(omega * mu0 / (beta * (d + 1/p + 1/q) * (h**2 + p**2)))
        x = -10.5 * d
        while x < 9.5 * d:
            XX.append(x)
            if x > 0:
                Prf.append(c_amp * math.exp(-p * x))
            elif -d <= x <= 0:
                Prf.append(c_amp * (math.cos(h * x) - (p / h) * math.sin(h * x)))
            else:
                Prf.append(c_amp * (math.cos(h * d) + (p / h) * math.sin(h * d)) * math.exp(q * (x + d)))
            x += dx

    elif mode == "TM":
        pb    = (n2**2 / n1**2) * p
        qb    = (n2**2 / n3**2) * q
        d_eff = ((pb**2 + h**2) / pb**2
                 * (d / n2**2
                    + (p**2 + h**2) / (pb**2 + h**2) / (n1**2 * p)
                    + (q**2 + h**2) / (qb**2 + h**2) / (n3**2 * q)))
        c_amp = -2 * math.sqrt(omega * eps0 / (beta * d_eff))
        x = -10.5 * d
        while x < 9.5 * d:
            XX.append(x)
            if x > 0:
                Prf.append(-(h / pb) * c_amp * math.exp(-p * x) * (beta / (omega * eps0 * n1**2)))
            elif -d <= x <= 0:
                Prf.append(c_amp * (-(h / pb) * math.cos(h * x) + math.sin(h * x)) * (beta / (omega * eps0 * n2**2)))
            else:
                Prf.append(-c_amp * ((h / pb) * math.cos(h * d) + math.sin(h * d)) * math.exp(q * (x + d)) * (beta / (omega * eps0 * n3**2)))
            x += dx

    ni = len(Prf)

    # Power integration for normalisation
    power = sum(Prf[i] ** 2 * dx for i in range(1, ni - 1))
    power += (Prf[0] * Prf[-1]) * 0.5 * dx
    power = math.sqrt(abs(power))
    if power == 0:
        power = 1.0

    with open(filepath, "w") as fh:
        fh.write(f"{ni}   0.0\n")
        for i in range(ni):
            xv = XX[i] + 0.5 * d
            yv = Prf[i] / power
            fh.write(f"{xv:.5e}  {yv:.5e}  0.0\n")


# ─────────────────────────────────────────────────────────────
# Sellmeier refractive index calculation
# ─────────────────────────────────────────────────────────────

def calc_sellmeier(layer, wl_nm, temp=25.0):
    """
    Calculate refractive index using Sellmeier equation.
    layer: 0=core, 1=over, 2=buff
    wl_nm: wavelength in nm
    temp : temperature in deg C

    Coefficients stored as flat list [B1,L1, B2,L2, B3,L3] (pairs).

    TypeA: n² = 1 + Σ Bi·λ²/(λ²−Li²)   (Li is squared)   — VB frm_RIndex TypeA
    TypeB: n² = 1 + Σ Bi·λ²/(λ²−Li)    (Li used directly) — VB frm_RIndex TypeB
    Table: use Calc_IndxTable()
    """
    import globals as g

    s = g.sellmeier[layer]
    wl_um = wl_nm / 1000.0  # nm → μm

    typ = s.get("type", "B")

    if typ in ("A", "TypeA"):
        c = s["coeffs"]
        n2 = 1.0
        for k in range(3):
            B = c[k * 2]     if k * 2     < len(c) else 0.0
            L = c[k * 2 + 1] if k * 2 + 1 < len(c) else 0.0
            if L != 0:
                n2 += B * wl_um ** 2 / (wl_um ** 2 - L ** 2)   # L squared
        n = math.sqrt(max(n2, 1.0))

    elif typ in ("B", "TypeB"):
        c = s["coeffs"]
        n2 = 1.0
        for k in range(3):
            B = c[k * 2]     if k * 2     < len(c) else 0.0
            L = c[k * 2 + 1] if k * 2 + 1 < len(c) else 0.0
            if L != 0:
                n2 += B * wl_um ** 2 / (wl_um ** 2 - L)        # L NOT squared
        n = math.sqrt(max(n2, 1.0))

    elif typ == "Table":
        n = Calc_IndxTable(layer, wl_nm)

    else:
        c = s.get("coeffs", [1.44])
        n = c[0] if c else 1.44

    # Temperature correction
    tc = g.nr_Tcoeff[layer]
    n += tc * (temp - 25.0)

    return n


# ─────────────────────────────────────────────────────────────
# .sel file I/O
# ─────────────────────────────────────────────────────────────

def load_sel_file(filepath):
    """
    Parse a Sellmeier coefficient file (.sel).

    File format (VB frm_RIndex.frm Read_Material convention):
        Line 1 : material name
        Line 2 : TypeA / TypeB / Table
        Line 3 : N  (number of coefficient pairs / table rows)
        Lines 4..4+N-1 : two values per line  "a  L"
        Remaining (if present):
            TypeB : WL_min  WL_max  T_coeff
            Table : (blank)  (blank)  T_coeff
            TypeA : (nothing)

    Returns dict:
        {"name": str, "type": str,
         "coeffs": [B1,L1,B2,L2,B3,L3],   # Sellmeier pairs (flat, up to 3)
         "table": [(wl_nm, n), ...],         # populated for Table type
         "wl_min": float, "wl_max": float,
         "tcoeff": float}
    """
    result = {
        "name": "", "type": "B",
        "coeffs": [0.0] * 6,
        "table": [],
        "wl_min": 0.0, "wl_max": 9999.0,
        "tcoeff": 0.0,
    }
    try:
        with open(filepath, "r", encoding="utf-8", errors="replace") as fh:
            raw = fh.readlines()

        lines = [l.rstrip() for l in raw]
        idx = 0

        def nextline():
            nonlocal idx
            while idx < len(lines):
                v = idx
                idx += 1
                return lines[v]
            return ""

        result["name"] = nextline().strip()
        result["type"] = nextline().strip()
        N = int(nextline().strip() or "0")

        coeffs = [0.0] * 6
        table  = []
        for i in range(N):
            parts = nextline().split()
            if len(parts) >= 2:
                a = float(parts[0])
                L = float(parts[1])
                if result["type"] == "Table":
                    table.append((a, L))       # wl_nm, n
                else:
                    if i < 3:
                        coeffs[i * 2]     = a
                        coeffs[i * 2 + 1] = L
        result["coeffs"] = coeffs
        result["table"]  = table

        # Remaining non-blank lines → WL_min, WL_max, T_coeff
        remaining = []
        while idx < len(lines):
            v = lines[idx].strip(); idx += 1
            if v:
                remaining.append(v)

        if len(remaining) >= 3:
            result["wl_min"] = _try_float(remaining[0])
            result["wl_max"] = _try_float(remaining[1])
            result["tcoeff"] = _try_float(remaining[2])
        elif len(remaining) == 2:
            result["wl_min"] = _try_float(remaining[0])
            result["wl_max"] = _try_float(remaining[1])
        elif len(remaining) == 1:
            result["tcoeff"] = _try_float(remaining[0])

    except Exception:
        pass
    return result


def _try_float(s):
    try:
        return float(s.split()[0])
    except Exception:
        return 0.0


def scan_sel_files(data_dir):
    """
    Scan *data_dir* for *.sel files.
    Returns list of (display_name, filepath, data_dict) sorted by name.
    """
    import os, glob
    items = []
    for fp in sorted(glob.glob(os.path.join(data_dir, "*.sel"))):
        d = load_sel_file(fp)
        name = d["name"] or os.path.splitext(os.path.basename(fp))[0]
        items.append((name, fp, d))
    return items


def apply_sel_to_globals(layer, sel_data):
    """
    Store parsed .sel data into global arrays for the given layer (0/1/2).
    Updates g.sellmeier[layer], g.nr_Tcoeff[layer], g.nr_type[layer],
    and (for Table type) g.nr_table / g.nr_index.
    """
    import globals as g

    g.sellmeier[layer]["name"]   = sel_data["name"]
    g.sellmeier[layer]["type"]   = sel_data["type"]
    g.sellmeier[layer]["coeffs"] = list(sel_data["coeffs"])
    g.nr_Tcoeff[layer]           = sel_data["tcoeff"]
    g.nr_type[layer]             = sel_data["type"]
    g.nr_Wlim[layer][0]         = sel_data["wl_min"]
    g.nr_Wlim[layer][1]         = sel_data["wl_max"]

    if sel_data["type"] == "Table":
        rows = sel_data["table"]
        g.nr_table[layer][0] = [r[0] for r in rows]   # wavelengths
        g.nr_table[layer][1] = [r[1] for r in rows]   # indices
        g.nr_index[layer]    = len(rows)
