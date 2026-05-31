"""
AWG WDM Design Utility - Global Variables and Constants
Ported from Visual Basic 5.0/6.0 to Python
Original: ETRI, (C) 2000 by Jang-Uk Shin and ETRI
"""

import math

# Program information
VER_STRING = "AWG WDM Design Utility V2.01 (Python)"
COPYRIGHT = "(C) Copyright by Jang-Uk Shin and ETRI (2000) - Python port"

# Constants
PI = math.pi
YES = 1
NO_ = 0

# Flags
f_saved = YES
f_completed = NO_
f_bpmcoeff = NO_
f_fsrcalc = NO_

# Program global variables
nc = 0          # number of input channels
na = 0          # number of array waveguides
mphase = 0.0    # phase order
wl_c = 0.0      # center wavelength [m]
wl_d = 0.0      # channel spacing [m]
wl_s = 0.0      # start wavelength [m]
wl_e = 0.0      # end wavelength [m]
wg_w = 0.0      # waveguide width [m]
wg_h = 0.0      # waveguide height [m]
pol = "X"       # polarization
n_over = 0.0    # over cladding index
n_core = 0.0    # core index
n_buff = 0.0    # buffer/cladding index
n_eff0 = 0.0    # waveguide mode index
n_eff1 = 0.0    # core region slab effective index
n_eff2 = 0.0    # cladding region effective index
d_aw = 0.0      # distance between array waveguides [m]
d_io = 0.0      # distance between IO waveguides [m]
d_path = 0.0    # path difference of array waveguides [m]
fsr = 0.0       # free spectral range [m]
f_l = 0.0       # focal length [m]
angle_d = 0.0   # angle between array waveguides [rad]
angle_i = 0.0   # angle between IO ports [rad]
betta = 0.0     # waveguide mode propagation constant
offs_1 = 0.0    # in/out waveguide offset
offs_2_0 = 0.0  # AWG waveguide offset A
offs_2_1 = 0.0  # AWG waveguide offset B
offs_2_2 = 0.0  # AWG waveguide offset C
fa_awg = 0.0    # AWG side slab focal length adjustment
fa_io = 0.0     # IO side slab focal length adjustment
gap_sw = 0.0    # gap between slab and waveguides
cf_Gaus = 0.35  # Gaussian coupling coefficient width

l_total = 0.0   # total length of device
d_port = 0.0    # port distance
d_slab = 0.0    # distance between slabs
a_slab = 0.0    # slab angle
l_S_min = 0.0   # straight line minimum length
l_B_min = 0.0   # minimum bend radius
nDum_I = 0      # input dummy waveguides
nDum_O = 0      # output dummy waveguides

# Taper/parabolic waveguide parameters
nchk_in_t = 0
nchk_in_p = 0
W_in_item = 0.0
L_in_item = 0.0
nchk_in_s = 0          # Segment checkbox
N_in_seg  = 4          # N segment
Por_in_seg = 80.0      # Portion (%)

nchk_out_t = 0
nchk_out_p = 0
W_out_item = 0.0
L_out_item = 0.0
nchk_out_s = 0
N_out_seg  = 4
Por_out_seg = 80.0

nchk_Sin_t = 0
nchk_Sin_p = 0
W_Sin_item = 0.0
L_Sin_item = 0.0
nchk_Sin_s = 0
N_Sin_seg  = 4
Por_Sin_seg = 80.0

nchk_Sout_t = 0
nchk_Sout_p = 0
W_Sout_item = 0.0
L_Sout_item = 0.0
nchk_Sout_s = 0
N_Sout_seg  = 4
Por_Sout_seg = 80.0

# AWG geometry arrays (1-indexed, index 0 unused)
AX = []   # AWG waveguide start x
AY = []   # AWG waveguide start y
BX = []   # AWG waveguide end x
BY = []   # AWG waveguide end y
OX = []   # AWG arc center x
OY = []   # AWG arc center y
R = []    # AWG arc radius
AN = []   # AWG waveguide angles

DDX = []  # input slab port x
DDY = []  # input slab port y
EX = []   # input slab arc x
EY = []   # input slab arc y
FX = []   # output slab port x
FY = []   # output slab port y
GX = []   # output slab arc x
GY = []   # output slab arc y
cx = []   # output arc center x
CY = []   # output arc center y
AA = []   # output port angles

OF_AW = []      # AWG waveguide offset
Cpl_Coeff = []  # coupling coefficient amplitudes
CplpCoeff = []  # coupling coefficient phases

Awg_Loss = []    # [dB/cm]
Awg_Phase = []   # [degrees]
Awg_Length = []  # [um]
Awg_Arc_R = []   # [um]
Awg_Arc_L = []   # [um]
Awg_Width = []   # [um]
Awg_Height = []  # [um]

ncc = 0.0  # group index related parameter
chk_1N_value = True  # True = 1xN configuration (used by simulation)

# Refractive index tables: nr_table[layer][point][col]
# layer: 0=core, 1=over, 2=buff, 3=spare, 4=mode_table
# point: data point index (1-based in VB, stored 1..n here)
# col: 1=wavelength, 2=index
nr_table = [[[] for _ in range(3)] for _ in range(5)]  # [layer][col] = list of values
nr_index = [0] * 5   # number of entries per layer
nr_type = [""] * 5   # type: "A", "B", or "T" (table)
nr_Wlim = [[0.0, 0.0, 0.0] for _ in range(5)]  # wavelength limits
nr_Tcoeff = [0.0] * 5  # temperature coefficient

# Sellmeier coefficients for refractive index calculation
# Stored per material (0=core, 1=over, 2=buff)
# Type A: n^2 = A + B*wl^2/(wl^2-C) + D*wl^2/(wl^2-E)  (wl in um)
# Type B: n^2 = 1 + B1*wl^2/(wl^2-C1) + B2*wl^2/(wl^2-C2) + B3*wl^2/(wl^2-C3)
sellmeier = [
    {"name": "", "type": "B", "coeffs": [0.0]*6},  # core
    {"name": "", "type": "B", "coeffs": [0.0]*6},  # over
    {"name": "", "type": "B", "coeffs": [0.0]*6},  # buff
]

# Temperature for refractive index calculation
temperature = 25.0  # degrees C

def arcsin(x):
    """Arc sine function."""
    return math.atan(x / math.sqrt(-x * x + 1))

def Max(a, b):
    return a if a >= b else b

def Min(a, b):
    return a if a <= b else b

def init_arrays(new_nc, new_na):
    """Initialize/resize global arrays when nc or na changes."""
    global AX, AY, BX, BY, OX, OY, R, AN
    global DDX, DDY, EX, EY, FX, FY, GX, GY, cx, CY, AA
    global OF_AW, Cpl_Coeff, CplpCoeff
    global Awg_Loss, Awg_Phase, Awg_Length, Awg_Arc_R, Awg_Arc_L, Awg_Width, Awg_Height

    n = max(new_na, 1)
    AX = [0.0] * (n + 1)
    AY = [0.0] * (n + 1)
    BX = [0.0] * (n + 1)
    BY = [0.0] * (n + 1)
    OX = [0.0] * (n + 1)
    OY = [0.0] * (n + 1)
    R  = [0.0] * (n + 1)
    AN = [0.0] * (n + 1)
    OF_AW = [0.0] * (n + 1)
    Cpl_Coeff = [0.0] * (n + 1)
    CplpCoeff = [0.0] * (n + 1)
    Awg_Loss   = [0.0] * (n + 1)
    Awg_Phase  = [0.0] * (n + 1)
    Awg_Length = [0.0] * (n + 1)
    Awg_Arc_R  = [0.0] * (n + 1)
    Awg_Arc_L  = [0.0] * (n + 1)
    Awg_Width  = [0.0] * (n + 1)
    Awg_Height = [0.0] * (n + 1)

    nc2 = max(new_nc, 1)
    DDX = [0.0] * (nc2 + 1)
    DDY = [0.0] * (nc2 + 1)
    EX  = [0.0] * (nc2 + 1)
    EY  = [0.0] * (nc2 + 1)
    FX  = [0.0] * (nc2 + 1)
    FY  = [0.0] * (nc2 + 1)
    GX  = [0.0] * (nc2 + 1)
    GY  = [0.0] * (nc2 + 1)
    cx  = [0.0] * (nc2 + 1)
    CY  = [0.0] * (nc2 + 1)
    AA  = [0.0] * (nc2 + 1)
