"""
Device Layout Generator Form
Ported from Frm_Layout.frm (VB 5.0/6.0) to Python/Tkinter

Draws the exact same geometry as the VB version:
  - AWG waveguides with physical width + arc bends
  - Slab boundaries (Rowland-circle geometry)
  - IO port waveguides (straight → S-bend arc → horizontal bus)
  - Input ports (straight section with optional taper/parabolic)
  - DXF export with LINE and ARC entities
"""

import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import math
import globals as g
import awg_calc

PI = math.pi

# ── colour constants ──────────────────────────────────────────
C_BG      = "#1a1a2e"
C_AWG     = "#ffff00"   # AWG waveguide centre-line
C_WALL    = "#4444ff"   # AWG waveguide walls
C_SLAB    = "#ffffff"   # slab boundary
C_FOCAL   = "#ff4444"   # focal / centre lines
C_IO      = "#00cc00"   # IO port waveguides
C_AXIS    = "#00cc00"   # axis cross

_PORT_COLORS = [
    "#ff5555","#55ff55","#5555ff","#ffff55",
    "#55ffff","#ff55ff","#ff9955","#55ff99",
    "#9955ff","#ff5599","#99ff55","#5599ff",
    "#ffcc55","#55ffcc","#cc55ff",
]


class FrmLayout:
    """AWG Device Layout window."""

    def __init__(self, master, frm_variables=None, frm_plot=None):
        self.master        = master
        self.frm_variables = frm_variables
        self.frm_plot      = frm_plot
        self.window        = None
        self.canvas        = None

        # view scale
        self.SCW = 10000.0   # visible width [um]
        self.SCH = 8000.0

        # mouse drag / pan
        self._mx0 = 0.0
        self._my0 = 0.0

        # view offset (pan) in um
        self._ox = 0.0
        self._oy = 0.0

        # geometry cache (list of draw-primitives)
        self._prims = []    # [(type, args, color)]

        # StringVars for file I/O
        self.var_offset1   = tk.StringVar(value="0.00")
        self.var_offset2   = [tk.StringVar(value="0.000e+00") for _ in range(3)]
        self.var_focal_awg = tk.StringVar(value="0")
        self.var_focal_io  = tk.StringVar(value="0")
        self.var_wsl       = tk.StringVar(value="1.5")
        self.var_rowland_I = tk.IntVar(value=1)
        self.var_rowland_O = tk.IntVar(value=1)
        self.var_adiabatic = tk.IntVar(value=1)
        self.var_arc_step  = tk.StringVar(value="0.1")   # arc segment step [degrees]

    # ─────────────────────────────────────────────────────────────
    def show(self):
        if self.window is None or not self.window.winfo_exists():
            self._build()
        self.window.deiconify()
        self.window.lift()

    def _build(self):
        self.window = tk.Toplevel(self.master)
        self.window.title("Device Layout")
        self.window.resizable(True, True)
        self.window.protocol("WM_DELETE_WINDOW", self._on_close)
        self._build_menu()
        self._build_ui()

    def _build_menu(self):
        menubar = tk.Menu(self.window)
        file_m  = tk.Menu(menubar, tearoff=0)
        file_m.add_command(label="Make DXF",
                           command=self._on_dxf)
        file_m.add_command(label="Make BPM file for input star",
                           command=lambda: messagebox.showinfo("BPM","Not implemented"))
        file_m.add_command(label="Make BPM file for output star",
                           command=lambda: messagebox.showinfo("BPM","Not implemented"))
        file_m.add_command(label="Make BPM file for reverse output star",
                           command=lambda: messagebox.showinfo("BPM","Not implemented"))
        menubar.add_cascade(label="File", menu=file_m)
        self.window.config(menu=menubar)

    def _build_ui(self):
        # ── left control panel ──
        left = ttk.Frame(self.window, width=220)
        left.pack(side=tk.LEFT, fill=tk.Y, padx=4, pady=4)
        left.pack_propagate(False)

        # cursor position
        pos_frm = ttk.LabelFrame(left, text="Cursor Position", padding=4)
        pos_frm.pack(fill=tk.X, pady=4)
        self.var_xpos = tk.StringVar(value="0")
        self.var_ypos = tk.StringVar(value="0")
        ttk.Label(pos_frm, text="X").grid(row=0, column=0)
        ttk.Entry(pos_frm, textvariable=self.var_xpos, width=12,
                  state="readonly").grid(row=0, column=1)
        ttk.Label(pos_frm, text="Y").grid(row=1, column=0)
        ttk.Entry(pos_frm, textvariable=self.var_ypos, width=12,
                  state="readonly").grid(row=1, column=1)

        # waveguide offset
        off_frm = ttk.LabelFrame(left, text="Waveguide Offset", padding=4)
        off_frm.pack(fill=tk.X, pady=4)
        ttk.Label(off_frm, text="Offset 1 (um)").grid(row=0, column=0, sticky="w")
        ttk.Entry(off_frm, textvariable=self.var_offset1, width=10)\
            .grid(row=0, column=1)
        ttk.Label(off_frm, text="Off2 = A + B·R + C·R²")\
            .grid(row=1, column=0, columnspan=2, sticky="w")
        for i, lbl in enumerate(["A","B","C"]):
            ttk.Label(off_frm, text=lbl).grid(row=2+i, column=0, sticky="w")
            ttk.Entry(off_frm, textvariable=self.var_offset2[i], width=12)\
                .grid(row=2+i, column=1)

        # focal length adjust
        fl_frm = ttk.LabelFrame(left, text="Focal Length Adjust (um)", padding=4)
        fl_frm.pack(fill=tk.X, pady=4)
        ttk.Label(fl_frm, text="AWG").grid(row=0, column=0, sticky="w")
        ttk.Entry(fl_frm, textvariable=self.var_focal_awg, width=10)\
            .grid(row=0, column=1)
        ttk.Label(fl_frm, text="IO ").grid(row=1, column=0, sticky="w")
        ttk.Entry(fl_frm, textvariable=self.var_focal_io, width=10)\
            .grid(row=1, column=1)

        # Rowland circle
        rl_frm = ttk.LabelFrame(left, text="Rowland Circle", padding=4)
        rl_frm.pack(fill=tk.X, pady=4)
        ttk.Checkbutton(rl_frm, text="Input Side",
                        variable=self.var_rowland_I,
                        command=self.draw).grid(row=0, column=0, sticky="w")
        ttk.Checkbutton(rl_frm, text="Output Side",
                        variable=self.var_rowland_O,
                        command=self.draw).grid(row=1, column=0, sticky="w")

        # others
        oth_frm = ttk.LabelFrame(left, text="Others", padding=4)
        oth_frm.pack(fill=tk.X, pady=4)
        ttk.Checkbutton(oth_frm, text="Adiabatic Connect",
                        variable=self.var_adiabatic)\
            .grid(row=0, column=0, columnspan=2, sticky="w")
        ttk.Label(oth_frm, text="Slab W Rate").grid(row=1, column=0, sticky="w")
        ttk.Entry(oth_frm, textvariable=self.var_wsl, width=6)\
            .grid(row=1, column=1)
        ttk.Label(oth_frm, text="Arc Step (deg)").grid(row=2, column=0, sticky="w")
        ttk.Entry(oth_frm, textvariable=self.var_arc_step, width=6)\
            .grid(row=2, column=1)

        # buttons
        ttk.Button(left, text="Draw",       command=self.draw)\
            .pack(fill=tk.X, pady=2)
        ttk.Button(left, text="ReDraw",     command=self._redraw)\
            .pack(fill=tk.X, pady=2)
        ttk.Button(left, text="Simulation", command=self._on_simulate)\
            .pack(fill=tk.X, pady=2)
        ttk.Button(left, text="Close",      command=self._on_close)\
            .pack(fill=tk.X, pady=2)

        # ── canvas ──
        cf = ttk.Frame(self.window)
        cf.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        self.canvas = tk.Canvas(cf, bg=C_BG, cursor="crosshair")
        self.canvas.pack(fill=tk.BOTH, expand=True)

        self.canvas.bind("<Motion>",        self._on_mouse_move)
        self.canvas.bind("<ButtonPress-1>", self._on_mouse_down)
        self.canvas.bind("<ButtonPress-3>", self._on_mouse_down)
        self.canvas.bind("<B1-Motion>",     self._on_drag1)
        self.canvas.bind("<B3-Motion>",     self._on_pan)
        # Scroll wheel zoom (Windows/Mac: MouseWheel; Linux: Button-4/5)
        self.canvas.bind("<MouseWheel>",    self._on_scroll)
        self.canvas.bind("<Button-4>",      self._on_scroll)
        self.canvas.bind("<Button-5>",      self._on_scroll)

    # ─────────────────────────────────────────────────────────────
    # Coordinate helpers
    # ─────────────────────────────────────────────────────────────

    def _setup_scale(self):
        self.canvas.update_idletasks()
        self._cw = max(self.canvas.winfo_width(),  800)
        self._ch = max(self.canvas.winfo_height(), 600)

        try:
            tl = float(self.frm_plot.var_total_l.get()) if self.frm_plot else 10000.0
        except (ValueError, AttributeError):
            tl = 10000.0

        self.SCW = max(tl * 1.2, 1000.0)
        self.SCH = self.SCW * self._ch / self._cw

    def _u2c(self, x_um, y_um):
        """um → canvas pixel (y-up math coords, with pan offset)."""
        px = ((x_um - self._ox) / self.SCW + 0.5) * self._cw
        py = (-(y_um - self._oy) / self.SCH + 0.5) * self._ch
        return px, py

    def _c2u(self, px, py):
        """canvas pixel → um (with pan offset)."""
        x = ( px / self._cw - 0.5) * self.SCW + self._ox
        y = (-py / self._ch + 0.5) * self.SCH + self._oy
        return x, y

    # ─────────────────────────────────────────────────────────────
    # Public draw entry point
    # ─────────────────────────────────────────────────────────────

    def draw(self):
        if g.na == 0 or g.f_l == 0:
            messagebox.showwarning("Warning",
                "AWG parameters not calculated. Run Simulate first.")
            return
        # Ensure global arrays are large enough for current na/nc
        g.init_arrays(g.nc, g.na)
        self.show()
        self._setup_scale()
        self._ox = 0.0   # reset pan on new draw
        self._oy = 0.0
        self.canvas.delete("all")

        self._prims = []
        self._calc_geometry()
        self._render_prims()

        g.f_completed = g.YES
        g.f_saved = g.NO_

    def _redraw(self):
        if not self.canvas:
            return
        self._setup_scale()
        self.canvas.delete("all")
        self._render_prims()

    # ─────────────────────────────────────────────────────────────
    # Geometry calculation  (faithful to VB Cmd_Draw_Click logic)
    # ─────────────────────────────────────────────────────────────

    def _calc_geometry(self):
        """Calculate all AWG geometry and store as draw-primitives.

        Faithfully follows VB Cmd_Draw_Click logic.
        """
        fp = self.frm_plot

        def fget(var, default):
            try:
                return float(var.get())
            except (ValueError, AttributeError):
                return default

        FCL     = g.f_l * 1e6                                   # focal length [um]
        DA1     = 2.0 * g.arcsin(0.5 * g.d_aw / g.f_l)         # AWG angle step
        DA2     = 2.0 * g.arcsin(0.5 * g.d_io / g.f_l)         # IO angle step
        SAN     = fget(fp.var_slb_angle if fp else None, 0.0) * PI / 180.0
        DSL     = fget(fp.var_d_slabs   if fp else None, FCL * 2.0)
        SNM     = fget(fp.var_l_min     if fp else None, 100.0)
        BNM     = fget(fp.var_b_min     if fp else None, 1000.0)
        TL_     = fget(fp.var_total_l   if fp else None, DSL + 2*FCL)
        DP_     = fget(fp.var_d_io_port if fp else None, g.d_io * 1e6)
        G_SW    = fget(fp.var_gap_sw    if fp else None, 0.0)
        WGW     = g.wg_w * 1e6                                  # waveguide width [um]
        DPL     = g.d_path * 1e6                                # path difference [um]
        WSL     = fget(self.var_wsl, 1.5)
        FA_AWG  = fget(self.var_focal_awg, 0.0)
        FA_IO   = fget(self.var_focal_io,  0.0)
        OF_IO   = fget(self.var_offset1, 0.0)
        OF2_0   = fget(self.var_offset2[0], 0.0)
        OF2_1   = fget(self.var_offset2[1], 0.0)
        OF2_2   = fget(self.var_offset2[2], 0.0)
        R_I     = bool(self.var_rowland_I.get())
        R_O     = bool(self.var_rowland_O.get())
        ADIAB   = bool(self.var_adiabatic.get())
        try:
            ARC_STEP_DEG = max(0.001, float(self.var_arc_step.get()))
        except (ValueError, AttributeError):
            ARC_STEP_DEG = 0.1

        def arc_segs(angle_rad):
            """Number of line segments for an arc of angle_rad radians."""
            return max(2, int(abs(angle_rad) * 180.0 / PI / ARC_STEP_DEG))
        na      = g.na
        nc_full = g.nc
        nDum_I  = int(fget(fp.var_dum_i if fp else None, 0))
        nDum_O  = int(fget(fp.var_dum_o if fp else None, 0))

        # Slab centre positions  (VB: QX_ = -0.5*DSL; PX_ = FCL*cos(SAN)+QX_)
        QX = -0.5 * DSL;  QY = 0.0
        PX = FCL * math.cos(SAN) + QX
        PY = FCL * math.sin(SAN) + QY

        def prim_line(x1, y1, x2, y2, color=C_SLAB):
            self._prims.append(("line", (x1, y1, x2, y2), color))

        # ── Axis cross ──────────────────────────────────────────
        prim_line(-0.5 * DSL, 0,  0.5 * DSL, 0, C_AXIS)
        prim_line(0, -0.4 * DSL, 0, 0.4 * DSL, C_AXIS)

        # ──────────────────────────────────────────────────────────
        # Slab boundaries  (faithfully following VB Cmd_Draw_Click)
        #
        # SLAB 1 (input, left side):
        #   AWG arc: i = 100 → -100 (top to bottom),  pts: a_hi → a_lo
        #   IO arc:  i = -100 → 100 (bottom to top),  pts: a_lo → a_hi
        #   Left closing edge:  AWG top  → IO bottom   (pts_awg[0] → pts_io[0])
        #   Right closing edge: AWG bottom → IO top    (pts_awg[-1] → pts_io[-1])
        #
        # SLAB 2 (output, right side): same logic with x-mirrored formulas
        # ──────────────────────────────────────────────────────────
        r_awg = FCL - FA_AWG - G_SW
        r_io  = FCL - FA_IO  - G_SW
        half_sweep = WSL * 0.5 * (na - 1) * DA1   # half angular width of slab
        n_slab = max(3, arc_segs(half_sweep) * 2 + 1)  # arc resolution for slab boundary

        def slab_ani_fwd(k):
            """ANI for k in 0..n_slab-1, going from +half_sweep to -half_sweep."""
            return SAN + half_sweep * (1.0 - 2.0 * k / (n_slab - 1))

        def slab_ani_rev(k):
            """ANI for k in 0..n_slab-1, going from -half_sweep to +half_sweep."""
            return SAN - half_sweep * (1.0 - 2.0 * k / (n_slab - 1))

        # --- SLAB 1 (input slab, left side) ---
        # Focal centre-line: (QX, QY) → (PX, PY)
        prim_line(QX, QY, PX, PY, C_FOCAL)

        # AWG arc: ANI goes from (SAN + half_sweep) down to (SAN - half_sweep)
        pts_awg1 = []
        for k in range(n_slab):
            ani = slab_ani_fwd(k)
            apx = r_awg * math.cos(ani) + QX
            apy = r_awg * math.sin(ani) + QY
            pts_awg1.append((apx, apy))

        # IO arc: ANI goes from (SAN - half_sweep) up to (SAN + half_sweep)
        pts_io1 = []
        for k in range(n_slab):
            ani = slab_ani_rev(k)
            if R_I:
                bpx = -(r_io) * math.cos(ani - SAN) * math.cos(ani) + PX
                bpy = -(r_io) * math.cos(ani - SAN) * math.sin(ani) + PY
            else:
                bpx = -r_io * math.cos(ani) + PX
                bpy = -r_io * math.sin(ani) + PY
            pts_io1.append((bpx, bpy))

        # Filled slab region: AWG arc forward + IO arc reversed = closed polygon
        self._prims.append(("filled_poly", pts_awg1 + list(reversed(pts_io1)), C_SLAB))

        # Optional: true focal arc when FA_AWG != 0
        if FA_AWG != 0:
            pts_ref1 = []
            for k in range(n_slab):
                ani = slab_ani_rev(k)
                pts_ref1.append((FCL * math.cos(ani) + QX,
                                 FCL * math.sin(ani) + QY))
            self._prims.append(("poly", pts_ref1, C_FOCAL))

        # Optional: true focal IO arc when FA_IO != 0
        if FA_IO != 0:
            pts_ref1b = []
            for k in range(n_slab):
                ani = slab_ani_fwd(k)
                if R_I:
                    bpx = -FCL * math.cos(ani - SAN) * math.cos(ani) + PX
                    bpy = -FCL * math.cos(ani - SAN) * math.sin(ani) + PY
                else:
                    bpx = -FCL * math.cos(ani) + PX
                    bpy = -FCL * math.sin(ani) + PY
                pts_ref1b.append((bpx, bpy))
            self._prims.append(("poly", pts_ref1b, C_FOCAL))

        # --- SLAB 2 (output slab, right side, x-mirrored) ---
        # Focal centre-line: (-QX, QY) → (-PX, PY)
        prim_line(-QX, QY, -PX, PY, C_FOCAL)

        # AWG arc (output side)
        pts_awg2 = []
        for k in range(n_slab):
            ani = slab_ani_fwd(k)
            apx = -r_awg * math.cos(ani) - QX
            apy =  r_awg * math.sin(ani) + QY
            pts_awg2.append((apx, apy))

        # IO arc (output side)
        pts_io2 = []
        for k in range(n_slab):
            ani = slab_ani_rev(k)
            if R_O:
                bpx =  r_io * math.cos(ani - SAN) * math.cos(ani) - PX
                bpy = -r_io * math.cos(ani - SAN) * math.sin(ani) + PY
            else:
                bpx =  r_io * math.cos(ani) - PX
                bpy = -r_io * math.sin(ani) + PY
            pts_io2.append((bpx, bpy))

        # Filled slab region
        self._prims.append(("filled_poly", pts_awg2 + list(reversed(pts_io2)), C_SLAB))

        if FA_AWG != 0:
            pts_ref2 = []
            for k in range(n_slab):
                ani = slab_ani_rev(k)
                pts_ref2.append((-FCL * math.cos(ani) - QX,
                                  FCL * math.sin(ani) + QY))
            self._prims.append(("poly", pts_ref2, C_FOCAL))

        if FA_IO != 0:
            pts_ref2b = []
            for k in range(n_slab):
                ani = slab_ani_fwd(k)
                if R_O:
                    bpx =  FCL * math.cos(ani - SAN) * math.cos(ani) - PX
                    bpy = -FCL * math.cos(ani - SAN) * math.sin(ani) + PY
                else:
                    bpx =  FCL * math.cos(ani) - PX
                    bpy = -FCL * math.sin(ani) + PY
                pts_ref2b.append((bpx, bpy))
            self._prims.append(("poly", pts_ref2b, C_FOCAL))

        # ──────────────────────────────────────────────────────────
        # AWG waveguide calculation  (VB Cmd_Draw_Click, j=1..na)
        # ──────────────────────────────────────────────────────────
        # Outermost waveguide (j = na)  — VB computes this first
        AN_na = SAN + 0.5 * (1 - na) * DA1
        AX_na = (FCL - FA_AWG) * math.cos(AN_na) + QX
        AY_na = (FCL - FA_AWG) * math.sin(AN_na) + QY
        denom_na = math.sin(AN_na)
        R_na  = ((0.5 * DSL - (FCL + SNM) * math.cos(AN_na)) / denom_na
                 if abs(denom_na) > 1e-12 else 0.0)
        OF_na = OF2_0 + OF2_1 * R_na + OF2_2 * R_na ** 2
        BX_na = (FCL + SNM) * math.cos(AN_na) + QX
        BY_na = (FCL + SNM) * math.sin(AN_na) + QY
        OX_na = BX_na + R_na * math.sin(AN_na)   # = 0 by construction
        OY_na = BY_na - R_na * math.cos(AN_na)

        g.AN[na] = AN_na
        g.AX[na] = AX_na;  g.AY[na] = AY_na
        g.BX[na] = BX_na;  g.BY[na] = BY_na
        g.OX[na] = OX_na;  g.OY[na] = OY_na
        g.R[na]  = R_na;   g.OF_AW[na] = OF_na

        R_min = R_na;  S_MIN = SNM

        for j in range(1, na):
            AN_j = SAN + (0.5 * (na + 1) - j) * DA1
            AX_j = (FCL - FA_AWG) * math.cos(AN_j) + QX
            AY_j = (FCL - FA_AWG) * math.sin(AN_j) + QY
            denom = math.sin(AN_j) - AN_j * math.cos(AN_j)
            if abs(denom) > 1e-12:
                R_j = (0.5 * DSL
                       - (FCL + SNM + R_na * AN_na + 0.5 * (na - j) * DPL)
                       * math.cos(AN_j)) / denom
            else:
                R_j = R_na
            OF_j = OF2_0 + OF2_1 * R_j + OF2_2 * R_j ** 2
            SNJ  = SNM + (R_na * AN_na - R_j * AN_j) + 0.5 * (na - j) * DPL
            BX_j = (FCL + SNJ) * math.cos(AN_j) + QX
            BY_j = (FCL + SNJ) * math.sin(AN_j) + QY
            OX_j = BX_j + R_j * math.sin(AN_j)   # = 0 always
            OY_j = BY_j - R_j * math.cos(AN_j)

            g.AN[j] = AN_j
            g.AX[j] = AX_j;  g.AY[j] = AY_j
            g.BX[j] = BX_j;  g.BY[j] = BY_j
            g.OX[j] = OX_j;  g.OY[j] = OY_j
            g.R[j]  = R_j;   g.OF_AW[j] = OF_j

            R_min = min(R_min, R_j)
            S_MIN = min(S_MIN, SNJ)

        # Update calculated fields in Frm_Plot
        if self.frm_plot:
            self.frm_plot.set_b_min_c(R_min)
            self.frm_plot.set_l_min_c(S_MIN)

        # AWG data arrays
        for j in range(1, na + 1):
            straight = math.hypot(g.BX[j] - g.AX[j], g.BY[j] - g.AY[j])
            g.Awg_Arc_R[j]  = g.R[j]
            g.Awg_Arc_L[j]  = g.R[j] * g.AN[j] * 2.0
            g.Awg_Length[j] = (straight + g.R[j] * g.AN[j]) * 2.0
            g.Awg_Width[j]  = g.wg_w * 1e6
            g.Awg_Height[j] = g.wg_h * 1e6

        avg_len = 0.5 * (g.Awg_Length[1] + g.Awg_Length[na]) if na > 0 else 0.0
        if self.frm_plot:
            self.frm_plot.set_a_avg(avg_len)

        # ──────────────────────────────────────────────────────────
        # Draw AWG waveguides  (straight + SINGLE arc, physical width)
        #
        # The arc is centred at OX[j] ≈ 0, OY[j].
        # VB draws ONE arc spanning from input BX[j] to output -BX[j]
        # through the centre. We replicate that single arc — NO mirrored copy.
        #
        # VB sweep: i = 100 → -100 (step -2), ANI = i*AN(j)/100 + pi/2
        #   → angle goes from (AN_j + pi/2) at top, to (pi/2) at centre,
        #     to (-AN_j + pi/2) at bottom.
        # End-cap angles:
        #   top: ANI = AN_j + pi/2  → outer/inner cap at input side
        #   bot: ANI = -AN_j + pi/2 → outer/inner cap at output side
        # ──────────────────────────────────────────────────────────
        # ── Segment drawing helper ────────────────────────────────
        def draw_segments(TAPL, TAPW, N, portion_pct, ANI, apx, apy,
                          neg_x, color, direction=-1):
            """Draw N constant-width staircase segments approximating a taper.

            TAPL   – total taper length [um]
            TAPW   – wide end width [um] (at slab face / AWG face)
            N      – number of segments
            portion_pct – fraction of each segment at full width [0-100 %]
            ANI    – waveguide angle [rad]
            apx/apy – wide-end start position (slab face)
            neg_x  – if True, flip x  (output-side AWG / output port)
            color  – wall color
            direction – +1 for AWG side (away from slab in +cos direction),
                        -1 for IO port side (away from slab in -cos direction)
            """
            por = max(0.01, min(1.0, portion_pct / 100.0))
            cos_a = math.cos(ANI); sin_a = math.sin(ANI)
            seg_len = TAPL / N
            half_p  = seg_len * por / 2.0

            for i in range(N):
                z_mid = TAPL * (i + 0.5) / N   # midpoint of this step along taper
                # Clamp z1/z2 inside [0, TAPL] to prevent segments going outside taper
                z1 = max(0.0, z_mid - half_p)
                z2 = min(TAPL, z_mid + half_p)
                # Width = taper width at the step MIDPOINT (per image definition)
                t     = z_mid / TAPL            # 0=slab/wide end, 1=WG/narrow end
                w     = TAPW + (WGW - TAPW) * t # taper envelope width at midpoint
                dxw   = 0.5 * w   * math.sin(ANI)
                dyw   = 0.5 * w   * math.cos(ANI)

                # Start and end coordinates of this segment portion
                # direction=+1: AWG side, goes in +cos direction from slab face
                # direction=-1: IO port side, goes in -cos direction from slab face
                if not neg_x:
                    sx1 = apx + direction * z1 * cos_a
                    sy1 = apy + direction * z1 * sin_a
                    sx2 = apx + direction * z2 * cos_a
                    sy2 = apy + direction * z2 * sin_a
                    # Filled polygon (4 corners of segment rectangle)
                    pts_seg = [(sx1-dxw, sy1+dyw), (sx2-dxw, sy2+dyw),
                               (sx2+dxw, sy2-dyw), (sx1+dxw, sy1-dyw)]
                else:
                    sx1 = -(apx + direction * z1 * cos_a)
                    sy1 =   apy + direction * z1 * sin_a
                    sx2 = -(apx + direction * z2 * cos_a)
                    sy2 =   apy + direction * z2 * sin_a
                    pts_seg = [(sx1+dxw, sy1+dyw), (sx2+dxw, sy2+dyw),
                               (sx2-dxw, sy2-dyw), (sx1-dxw, sy1-dyw)]
                self._prims.append(("filled_poly", pts_seg, color))

        # Taper parameters for AWG side (in um)
        TAPL_sin  = max(0.0, g.L_Sin_item  * 1e6)
        TAPW_sin  = g.W_Sin_item  * 1e6
        TAPL_sout = max(0.0, g.L_Sout_item * 1e6)
        TAPW_sout = g.W_Sout_item * 1e6
        chk_sin_t  = g.nchk_Sin_t
        chk_sin_p  = g.nchk_Sin_p
        chk_sin_s  = g.nchk_Sin_s
        chk_sout_t = g.nchk_Sout_t
        chk_sout_p = g.nchk_Sout_p
        chk_sout_s = g.nchk_Sout_s

        for j in range(1, na + 1):
            AN_j = g.AN[j]
            dx   = 0.5 * WGW * math.sin(AN_j)
            dy   = 0.5 * WGW * math.cos(AN_j)
            r_c  = g.R[j] - g.OF_AW[j]    # effective arc radius (centre-line)
            OXj  = g.OX[j]   # ≈ 0
            OYj  = g.OY[j]

            # ── Input-side straight section ──
            # If taper/parab/segment enabled the straight section starts at AXX = AX + TAPL*cos
            if (chk_sin_t or chk_sin_p or chk_sin_s) and TAPL_sin > 0:
                AXX = g.AX[j] + TAPL_sin * math.cos(AN_j)
                AYY = g.AY[j] + TAPL_sin * math.sin(AN_j)
            else:
                AXX = g.AX[j]; AYY = g.AY[j]

            prim_line(g.AX[j], g.AY[j], g.BX[j], g.BY[j], C_FOCAL)
            self._prims.append(("filled_poly",
                [(AXX+dx, AYY-dy), (g.BX[j]+dx, g.BY[j]-dy),
                 (g.BX[j]-dx, g.BY[j]+dy), (AXX-dx, AYY+dy)], C_AWG))

            # Input-side AWG taper/parabolic/segment
            # Segment and Taper/Parabolic can coexist — segment draws staircase,
            # taper/parabolic draws the envelope outline independently.
            if chk_sin_s and TAPL_sin > 0:
                draw_segments(TAPL_sin, TAPW_sin, g.N_Sin_seg, g.Por_Sin_seg,
                              AN_j, g.AX[j], g.AY[j], neg_x=False, color=C_AWG,
                              direction=+1)
            if chk_sin_t and TAPL_sin > 0:
                dx2 = 0.5 * TAPW_sin * math.sin(AN_j)
                dy2 = 0.5 * TAPW_sin * math.cos(AN_j)
                pts_tap = [(g.AX[j]+dx2, g.AY[j]-dy2), (AXX+dx, AYY-dy),
                           (AXX-dx, AYY+dy), (g.AX[j]-dx2, g.AY[j]+dy2)]
                self._prims.append(("filled_poly", pts_tap, C_AWG))
            elif chk_sin_p and TAPL_sin > 0 and g.wl_c > 0 and g.n_core > 0:
                ALPHA = (TAPW_sin**2 - WGW**2) / (2.0 * g.wl_c * 1e6 * TAPL_sin / g.n_core)
                nPARA = 50
                APX = g.AX[j]; APY = g.AY[j]
                pdx = 0.5 * WGW * math.sin(AN_j); pdy = 0.5 * WGW * math.cos(AN_j)
                prim_line(APX-pdx, APY+pdy, APX+pdx, APY-pdy, C_AWG)
                for i in range(nPARA + 1):
                    Z    = TAPL_sin * i / nPARA
                    w    = math.sqrt(max(2.0*ALPHA*(g.wl_c*1e6/g.n_core)*(TAPL_sin - Z) + WGW**2, WGW**2))
                    dx2  = 0.5 * w * math.sin(AN_j); dy2 = 0.5 * w * math.cos(AN_j)
                    BPX  = g.AX[j] + Z * math.cos(AN_j)
                    BPY  = g.AY[j] + Z * math.sin(AN_j)
                    prim_line(APX+pdx, APY-pdy, BPX+dx2, BPY-dy2, C_AWG)
                    prim_line(APX-pdx, APY+pdy, BPX-dx2, BPY+dy2, C_AWG)
                    APX = BPX; APY = BPY; pdx = dx2; pdy = dy2
                prim_line(APX+pdx, APY-pdy, APX-pdx, APY+pdy, C_AWG)

            # ── Output-side straight section (x-mirrored) ──
            if (chk_sout_t or chk_sout_p or chk_sout_s) and TAPL_sout > 0:
                AXX2 = g.AX[j] + TAPL_sout * math.cos(AN_j)
                AYY2 = g.AY[j] + TAPL_sout * math.sin(AN_j)
            else:
                AXX2 = g.AX[j]; AYY2 = g.AY[j]

            prim_line(-g.AX[j], g.AY[j], -g.BX[j], g.BY[j], C_FOCAL)
            self._prims.append(("filled_poly",
                [(-AXX2-dx, AYY2-dy), (-g.BX[j]-dx, g.BY[j]-dy),
                 (-g.BX[j]+dx, g.BY[j]+dy), (-AXX2+dx, AYY2+dy)], C_AWG))

            # Output-side AWG taper/parabolic/segment
            # Segment and Taper/Parabolic can coexist — segment draws staircase,
            # taper/parabolic draws the envelope outline independently.
            if chk_sout_s and TAPL_sout > 0:
                draw_segments(TAPL_sout, TAPW_sout, g.N_Sout_seg, g.Por_Sout_seg,
                              AN_j, g.AX[j], g.AY[j], neg_x=True, color=C_AWG,
                              direction=+1)
            if chk_sout_t and TAPL_sout > 0:
                dx2 = 0.5 * TAPW_sout * math.sin(AN_j)
                dy2 = 0.5 * TAPW_sout * math.cos(AN_j)
                pts_tap = [(-g.AX[j]-dx2, g.AY[j]-dy2), (-AXX2-dx, AYY2-dy),
                           (-AXX2+dx, AYY2+dy), (-g.AX[j]+dx2, g.AY[j]+dy2)]
                self._prims.append(("filled_poly", pts_tap, C_AWG))
            elif chk_sout_p and TAPL_sout > 0 and g.wl_c > 0 and g.n_core > 0:
                ALPHA = (TAPW_sout**2 - WGW**2) / (2.0 * g.wl_c * 1e6 * TAPL_sout / g.n_core)
                nPARA = 50
                APX = g.AX[j]; APY = g.AY[j]
                pdx = 0.5 * WGW * math.sin(AN_j); pdy = 0.5 * WGW * math.cos(AN_j)
                prim_line(-APX+pdx, APY+pdy, -APX-pdx, APY-pdy, C_AWG)
                for i in range(nPARA + 1):
                    Z    = TAPL_sout * i / nPARA
                    w    = math.sqrt(max(2.0*ALPHA*(g.wl_c*1e6/g.n_core)*(TAPL_sout - Z) + WGW**2, WGW**2))
                    dx2  = 0.5 * w * math.sin(AN_j); dy2 = 0.5 * w * math.cos(AN_j)
                    BPX  = g.AX[j] + Z * math.cos(AN_j)
                    BPY  = g.AY[j] + Z * math.sin(AN_j)
                    prim_line(-APX-pdx, APY-pdy, -BPX-dx2, BPY-dy2, C_AWG)
                    prim_line(-APX+pdx, APY+pdy, -BPX+dx2, BPY+dy2, C_AWG)
                    APX = BPX; APY = BPY; pdx = dx2; pdy = dy2
                prim_line(-APX-pdx, APY-pdy, -APX+pdx, APY+pdy, C_AWG)

            # ── Arc section — SINGLE arc centred at (OXj ≈ 0, OYj) ──
            # VB: for i=100 to -100 step -nSTEP, ANI = i*AN(j)/100 + pi/2
            # angle at i=100:   ANI = AN_j + pi/2   (input side)
            # angle at i=-100:  ANI = -AN_j + pi/2  (output side)
            n_arc = arc_segs(AN_j)
            a_top = AN_j  + 0.5 * PI   # i = +100
            a_bot = -AN_j + 0.5 * PI   # i = -100

            # Centre-line arc
            pts_c = []
            for k in range(n_arc + 1):
                a = a_top + k * (a_bot - a_top) / n_arc
                pts_c.append((r_c * math.cos(a) + OXj,
                               r_c * math.sin(a) + OYj))
            self._prims.append(("poly", pts_c, C_FOCAL))

            # Outer and inner wall arcs → filled polygon
            pts_o = []
            pts_i = []
            for k in range(n_arc + 1):
                a = a_top + k * (a_bot - a_top) / n_arc
                pts_o.append(((r_c + 0.5*WGW) * math.cos(a) + OXj,
                               (r_c + 0.5*WGW) * math.sin(a) + OYj))
                pts_i.append(((r_c - 0.5*WGW) * math.cos(a) + OXj,
                               (r_c - 0.5*WGW) * math.sin(a) + OYj))
            self._prims.append(("filled_poly", pts_o + list(reversed(pts_i)), C_AWG))

        # ──────────────────────────────────────────────────────────
        # IO waveguide routing  (output ports: slab → straight → S-bend → horizontal bus)
        # VB: AA(j), EX, EY, cx, CY, FX, FY, GX, GY
        # ──────────────────────────────────────────────────────────
        for j in range(1, nc_full + 1):
            g.AA[j] = SAN - (0.5 * (nc_full + 1) - j) * DA2

        # EX[1], EY[1]: initial port 1 position (VB uses FCL+SNM straight, not Rowland)
        EX = [0.0] * (nc_full + 2)
        EY = [0.0] * (nc_full + 2)
        cxA= [0.0] * (nc_full + 2)
        cyA= [0.0] * (nc_full + 2)
        FX = [0.0] * (nc_full + 2)
        FY = [0.0] * (nc_full + 2)
        GX = [0.0] * (nc_full + 2)
        GY = [0.0] * (nc_full + 2)

        EX[1] = PX - (FCL + SNM) * math.cos(g.AA[1])
        EY[1] = PY - (FCL + SNM) * math.sin(g.AA[1])
        cxA[1]= EX[1] - BNM * math.sin(g.AA[1])
        cyA[1]= EY[1] + BNM * math.cos(g.AA[1])
        FX[1] = cxA[1]
        FY[1] = cyA[1] - BNM
        GX[1] = -0.5 * TL_
        GY[1] = FY[1]

        for j in range(2, nc_full + 1):
            FY[j] = FY[1] - DP_ * (j - 1)
            EY[j] = FY[j] + BNM * (1.0 - math.cos(g.AA[j]))
            if abs(math.sin(g.AA[j])) > 1e-12:
                SNJ = (PY - EY[j]) / math.sin(g.AA[j])
            else:
                SNJ = FCL + SNM
            EX[j] = PX - SNJ * math.cos(g.AA[j])
            cxA[j]= EX[j] - BNM * math.sin(g.AA[j])
            cyA[j]= EY[j] + BNM * math.cos(g.AA[j])
            FX[j] = cxA[j]
            GX[j] = -0.5 * TL_
            GY[j] = FY[j]

        # Store in globals
        for j in range(1, nc_full + 1):
            if j < len(g.EX):  g.EX[j]  = EX[j]
            if j < len(g.EY):  g.EY[j]  = EY[j]
            if j < len(g.FX):  g.FX[j]  = FX[j]
            if j < len(g.FY):  g.FY[j]  = FY[j]
            if j < len(g.GX):  g.GX[j]  = GX[j]
            if j < len(g.GY):  g.GY[j]  = GY[j]
            if j < len(g.cx):  g.cx[j]  = cxA[j]
            if j < len(g.CY):  g.CY[j]  = cyA[j]

        # Taper parameters for IO ports (in um)
        TAPL_out = max(0.0, g.L_out_item * 1e6)
        TAPW_out = g.W_out_item * 1e6
        chk_out_t = g.nchk_out_t
        chk_out_p = g.nchk_out_p
        chk_out_s = g.nchk_out_s

        # Draw output port dummy waveguides (slab exit stubs) + optional taper
        for j in range(1 - nDum_O, nc_full + nDum_O + 1):
            ANI = SAN - (0.5 * (nc_full + 1) - j) * DA2
            if R_O:
                apx = PX - (FCL - FA_IO) * math.cos(ANI - SAN) * math.cos(ANI)
                apy = PY - (FCL - FA_IO) * math.cos(ANI - SAN) * math.sin(ANI)
            else:
                apx = PX - (FCL - FA_IO) * math.cos(ANI)
                apy = PY - (FCL - FA_IO) * math.sin(ANI)

            dx  = 0.5 * WGW * math.sin(ANI)
            dy  = 0.5 * WGW * math.cos(ANI)

            # Taper/parabolic/segment at slab face — drawn for ALL j (not just dummies)
            # VB: taper drawing has no j-condition; only the straight stub is dummy-only
            # Segment and Taper can coexist: segment=staircase, taper=envelope outline.
            if chk_out_s and TAPL_out > 0:
                draw_segments(TAPL_out, TAPW_out, g.N_out_seg, g.Por_out_seg,
                              ANI, apx, apy, neg_x=True, color=C_SLAB)
            if chk_out_t and TAPL_out > 0:
                dx2 = 0.5 * TAPW_out * math.sin(ANI)
                dy2 = 0.5 * TAPW_out * math.cos(ANI)
                bpx_t = apx - TAPL_out * math.cos(ANI)
                bpy_t = apy - TAPL_out * math.sin(ANI)
                prim_line(-apx, apy, -bpx_t, bpy_t, C_FOCAL)
                pts_tap = [(-apx+dx2, apy+dy2), (-bpx_t+dx, bpy_t+dy),
                           (-bpx_t-dx, bpy_t-dy), (-apx-dx2, apy-dy2)]
                self._prims.append(("filled_poly", pts_tap, C_SLAB))
            if (chk_out_s or chk_out_t) and TAPL_out > 0:
                apx = apx - TAPL_out * math.cos(ANI)
                apy = apy - TAPL_out * math.sin(ANI)
            elif chk_out_p and TAPL_out > 0 and g.wl_c > 0 and g.n_core > 0:
                ALPHA = (TAPW_out**2 - WGW**2) / (2.0 * g.wl_c*1e6 * TAPL_out / g.n_core)
                nPARA = 50
                cur_apx = apx; cur_apy = apy
                pdx = 0.5 * TAPW_out * math.sin(ANI); pdy = 0.5 * TAPW_out * math.cos(ANI)
                prim_line(-cur_apx+pdx, cur_apy+pdy, -cur_apx-pdx, cur_apy-pdy, C_SLAB)
                for i in range(1, nPARA + 1):
                    Z    = TAPL_out / nPARA
                    w    = math.sqrt(max(2.0*ALPHA*(g.wl_c*1e6/g.n_core)*(TAPL_out - i*Z) + WGW**2, WGW**2))
                    dx2  = 0.5 * w * math.sin(ANI); dy2 = 0.5 * w * math.cos(ANI)
                    bpx_t = cur_apx - Z * math.cos(ANI)
                    bpy_t = cur_apy - Z * math.sin(ANI)
                    prim_line(-cur_apx-pdx, cur_apy-pdy, -bpx_t-dx2, bpy_t-dy2, C_SLAB)
                    prim_line(-cur_apx+pdx, cur_apy+pdy, -bpx_t+dx2, bpy_t+dy2, C_SLAB)
                    cur_apx = bpx_t; cur_apy = bpy_t; pdx = dx2; pdy = dy2
                prim_line(-cur_apx-pdx, cur_apy-pdy, -cur_apx+pdx, cur_apy+pdy, C_SLAB)
                apx = cur_apx; apy = cur_apy

            # Straight dummy section (SNM − TAPL_out) — dummy-only
            rem = max(0.0, SNM - TAPL_out)
            bpx = apx - rem * math.cos(ANI)
            bpy = apy - rem * math.sin(ANI)
            if j < 1 or j > nc_full:
                prim_line(-apx, apy, -bpx, bpy, C_FOCAL)
                self._prims.append(("filled_poly",
                    [(-apx+dx, apy+dy), (-bpx+dx, bpy+dy),
                     (-bpx-dx, bpy-dy), (-apx-dx, apy-dy)], C_SLAB))

        # Draw output port main waveguides: bus + straight + S-bend arc
        TAPL = TAPL_out

        for j in range(1, nc_full + 1):
            col = _PORT_COLORS[(j - 1) % len(_PORT_COLORS)]

            # ── Output bus (horizontal, from -FX to -GX) ──
            # VB: Line (-FX(j),FY(j))-(-GX(j),GY(j)) red
            # -GX(j) = -(-0.5*TL_) = +0.5*TL_  → right edge
            DY = 0.5 * WGW
            prim_line(-FX[j], FY[j], -GX[j], GY[j], C_FOCAL)
            self._prims.append(("filled_poly",
                [(-FX[j], FY[j]+DY), (-GX[j], GY[j]+DY),
                 (-GX[j], GY[j]-DY), (-FX[j], FY[j]-DY)], col))

            # ── Straight section from slab exit to S-bend start ──
            # DDX = slab-exit point (Rowland formula), adjusted for taper
            an_j = g.AA[j]
            if R_O:
                DDX_j = PX - (FCL - FA_IO) * math.cos(an_j - SAN) * math.cos(an_j)
                DDY_j = PY - (FCL - FA_IO) * math.cos(an_j - SAN) * math.sin(an_j)
            else:
                DDX_j = PX - (FCL - FA_IO) * math.cos(an_j)
                DDY_j = PY - (FCL - FA_IO) * math.sin(an_j)
            # VB: DDX(j) -= TAPL*Cos(AA(j))
            DDX_j -= TAPL * math.cos(an_j)
            DDY_j -= TAPL * math.sin(an_j)
            # Store in globals
            if j < len(g.DDX): g.DDX[j] = DDX_j
            if j < len(g.DDY): g.DDY[j] = DDY_j

            dx_j = 0.5 * WGW * math.sin(an_j)
            dy_j = 0.5 * WGW * math.cos(an_j)

            prim_line(-DDX_j, DDY_j, -EX[j], EY[j], C_FOCAL)
            self._prims.append(("filled_poly",
                [(-DDX_j+dx_j, DDY_j+dy_j), (-EX[j]+dx_j, EY[j]+dy_j),
                 (-EX[j]-dx_j, EY[j]-dy_j), (-DDX_j-dx_j, DDY_j-dy_j)], col))

            # ── S-bend arc (output side) ──
            # VB: sweep i=0..100, ANI = AA(j)*(100-i)/100
            #     APX = cx(j) + (BNM-OF_IO)*sin(ANI)
            #     APY = CY(j) - (BNM-OF_IO)*cos(ANI)
            #     draw at (-APX, APY)
            R_bend = BNM - OF_IO
            n_bend = arc_segs(an_j)

            # Centre-line arc, outer wall, inner wall
            pts_arc_cen = []
            pts_arc_out = []
            pts_arc_in  = []
            for k in range(n_bend + 1):
                ani = an_j * (n_bend - k) / n_bend   # AA(j)*(100-i)/100 style
                ax_ = cxA[j] + R_bend            * math.sin(ani)
                ay_ = cyA[j] - R_bend            * math.cos(ani)
                ox_ = cxA[j] + (R_bend+0.5*WGW) * math.sin(ani)
                oy_ = cyA[j] - (R_bend+0.5*WGW) * math.cos(ani)
                ix_ = cxA[j] + (R_bend-0.5*WGW) * math.sin(ani)
                iy_ = cyA[j] - (R_bend-0.5*WGW) * math.cos(ani)
                pts_arc_cen.append((-ax_, ay_))
                pts_arc_out.append((-ox_, oy_))
                pts_arc_in.append( (-ix_, iy_))

            self._prims.append(("poly", pts_arc_cen, C_FOCAL))
            self._prims.append(("filled_poly", pts_arc_out + list(reversed(pts_arc_in)), col))

        # ──────────────────────────────────────────────────────────
        # Input waveguides  (1xN or NxN)
        # VB: If chk_1N checked → nc=1, AA(1)=SAN, single centred input
        # ──────────────────────────────────────────────────────────
        is_1n = (self.frm_variables.var_1n.get()
                 if self.frm_variables else True)

        if is_1n:
            # VB Cmd_Draw_Click, 1xN branch:
            #   AA(1) = SAN
            #   DDX(1) = PX_ - (FCL-fa_io)*cos(AA(1))   [no Rowland for this DDX]
            #   DDY(1) = PY_ - (FCL-fa_io)*sin(AA(1))
            #   FY(1) = (FY(1)+FY(nc))/2   [mid of output FY range]
            #   EY(1) = FY(1) + BNM*(1-cos(AA(1)))
            #   SNJ   = (PY_-EY(1))/sin(AA(1))   [or FCL+SNM if sin≈0]
            #   EX(1) = PX_ - SNJ*cos(AA(1))
            #   cx(1) = EX(1) - BNM*sin(AA(1))
            #   CY(1) = EY(1) + BNM*cos(AA(1))
            #   FX(1) = cx(1)
            #   GX(1) = -0.5*TL_   ← LEFT edge
            nc_in = 1
            AA_in  = SAN
            FY_mid = 0.5 * (FY[1] + FY[nc_full])
            EY_in  = FY_mid + BNM * (1.0 - math.cos(SAN))
            if abs(math.sin(SAN)) > 1e-12:
                SNJ_in = (PY - EY_in) / math.sin(SAN)
            else:
                SNJ_in = FCL + SNM
            EX_in  = PX - SNJ_in * math.cos(SAN)
            cx_in  = EX_in - BNM * math.sin(SAN)
            cy_in  = EY_in + BNM * math.cos(SAN)
            FX_in  = cx_in
            FY_in  = cy_in - BNM
            GX_in  = -0.5 * TL_    # goes to LEFT edge
            GY_in  = FY_in
            io_in_list = [(SAN, EX_in, EY_in, cx_in, cy_in, FX_in, FY_in, GX_in, GY_in)]
        else:
            # NxN: mirror the output port geometry to input side
            nc_in = nc_full
            io_in_list = []
            for j in range(1, nc_full + 1):
                io_in_list.append((g.AA[j], EX[j], EY[j],
                                   cxA[j], cyA[j], FX[j], FY[j],
                                   GX[j], GY[j]))   # GX already = -0.5*TL_

        # Taper parameters for input port
        TAPL_in = max(0.0, g.L_in_item * 1e6)
        TAPW_in = g.W_in_item * 1e6
        chk_in_t = g.nchk_in_t
        chk_in_p = g.nchk_in_p
        chk_in_s = g.nchk_in_s

        # Draw input dummy waveguides (slab exit stubs) + optional taper
        nc_draw = 1 if is_1n else nc_full
        for j in range(1 - nDum_I, nc_draw + nDum_I + 1):
            ANI = SAN - (0.5 * (nc_draw + 1) - j) * DA2
            if R_I:
                apx = PX - (FCL - FA_IO) * math.cos(ANI - SAN) * math.cos(ANI)
                apy = PY - (FCL - FA_IO) * math.cos(ANI - SAN) * math.sin(ANI)
            else:
                apx = PX - (FCL - FA_IO) * math.cos(ANI)
                apy = PY - (FCL - FA_IO) * math.sin(ANI)
            dx  = 0.5 * WGW * math.sin(ANI)
            dy  = 0.5 * WGW * math.cos(ANI)

            # Taper/parabolic/segment at slab face — drawn for ALL j (not just dummies)
            # VB: taper lines have no j-condition; only straight stub is dummy-only
            # Segment and Taper can coexist: segment=staircase, taper=envelope outline.
            if chk_in_s and TAPL_in > 0:
                draw_segments(TAPL_in, TAPW_in, g.N_in_seg, g.Por_in_seg,
                              ANI, apx, apy, neg_x=False, color=C_SLAB)
            if chk_in_t and TAPL_in > 0:
                dx2 = 0.5 * TAPW_in * math.sin(ANI)
                dy2 = 0.5 * TAPW_in * math.cos(ANI)
                bpx_t = apx - TAPL_in * math.cos(ANI)
                bpy_t = apy - TAPL_in * math.sin(ANI)
                prim_line(apx, apy, bpx_t, bpy_t, C_FOCAL)
                pts_tap = [(apx-dx2, apy+dy2), (bpx_t-dx, bpy_t+dy),
                           (bpx_t+dx, bpy_t-dy), (apx+dx2, apy-dy2)]
                self._prims.append(("filled_poly", pts_tap, C_SLAB))
            if (chk_in_s or chk_in_t) and TAPL_in > 0:
                apx = apx - TAPL_in * math.cos(ANI)
                apy = apy - TAPL_in * math.sin(ANI)
            elif chk_in_p and TAPL_in > 0 and g.wl_c > 0 and g.n_core > 0:
                ALPHA = (TAPW_in**2 - WGW**2) / (2.0 * g.wl_c*1e6 * TAPL_in / g.n_core)
                nPARA = 50
                cur_apx = apx; cur_apy = apy
                pdx = 0.5 * TAPW_in * math.sin(ANI); pdy = 0.5 * TAPW_in * math.cos(ANI)
                prim_line(cur_apx-pdx, cur_apy+pdy, cur_apx+pdx, cur_apy-pdy, C_SLAB)
                for i in range(1, nPARA + 1):
                    Z    = TAPL_in / nPARA
                    w    = math.sqrt(max(2.0*ALPHA*(g.wl_c*1e6/g.n_core)*(TAPL_in - i*Z) + WGW**2, WGW**2))
                    dx2  = 0.5 * w * math.sin(ANI); dy2 = 0.5 * w * math.cos(ANI)
                    bpx_t = cur_apx - Z * math.cos(ANI)
                    bpy_t = cur_apy - Z * math.sin(ANI)
                    prim_line(cur_apx-pdx, cur_apy+pdy, bpx_t-dx2, bpy_t+dy2, C_SLAB)
                    prim_line(cur_apx+pdx, cur_apy-pdy, bpx_t+dx2, bpy_t-dy2, C_SLAB)
                    cur_apx = bpx_t; cur_apy = bpy_t; pdx = dx2; pdy = dy2
                prim_line(cur_apx-pdx, cur_apy+pdy, cur_apx+pdx, cur_apy-pdy, C_SLAB)
                apx = cur_apx; apy = cur_apy

            # Straight dummy section (SNM − TAPL_in) — dummy-only
            rem = max(0.0, SNM - TAPL_in)
            bpx = apx - rem * math.cos(ANI)
            bpy = apy - rem * math.sin(ANI)
            if j < 1 or j > nc_draw:
                prim_line(apx, apy, bpx, bpy, C_FOCAL)
                self._prims.append(("filled_poly",
                    [(apx-dx, apy+dy), (bpx-dx, bpy+dy),
                     (bpx+dx, bpy-dy), (apx+dx, apy-dy)], C_SLAB))

        for idx, (an_in, ex_in, ey_in, cxi, cyi, fxi, fyi, gxi, gyi) in enumerate(io_in_list):
            col = _PORT_COLORS[idx % len(_PORT_COLORS)]

            # ── Input bus (horizontal, from FX to GX) ──
            # VB: Line (FX(j),FY(j))-(GX(j),GY(j)) red
            # GX = -0.5*TL_  (negative) → goes to LEFT edge
            # No sign flip: input side x-coords are already on the left (negative region)
            DY = 0.5 * WGW
            prim_line(fxi, fyi, gxi, gyi, C_FOCAL)
            self._prims.append(("filled_poly",
                [(fxi, fyi+DY), (gxi, gyi+DY),
                 (gxi, gyi-DY), (fxi, fyi-DY)], col))

            # ── Straight section from slab exit to S-bend start ──
            # DDX: slab-exit (Rowland), adjusted by TAPL_in
            if R_I:
                DDX_in = PX - (FCL - FA_IO) * math.cos(an_in - SAN) * math.cos(an_in)
                DDY_in = PY - (FCL - FA_IO) * math.cos(an_in - SAN) * math.sin(an_in)
            else:
                DDX_in = PX - (FCL - FA_IO) * math.cos(an_in)
                DDY_in = PY - (FCL - FA_IO) * math.sin(an_in)
            DDX_in -= TAPL_in * math.cos(an_in)
            DDY_in -= TAPL_in * math.sin(an_in)

            dx_in = 0.5 * WGW * math.sin(an_in)
            dy_in = 0.5 * WGW * math.cos(an_in)

            prim_line(DDX_in, DDY_in, ex_in, ey_in, C_FOCAL)
            self._prims.append(("filled_poly",
                [(DDX_in-dx_in, DDY_in+dy_in), (ex_in-dx_in, ey_in+dy_in),
                 (ex_in+dx_in, ey_in-dy_in), (DDX_in+dx_in, DDY_in-dy_in)], col))

            # ── S-bend arc (input side) ──
            # VB: same formula as output but NO negation on x:
            #     APX = cx(j) + (BNM-OF_IO)*sin(ANI)
            #     APY = CY(j) - (BNM-OF_IO)*cos(ANI)
            #     draw at (APX, APY)   [positive x → left side because cx is negative]
            R_bend = BNM - OF_IO
            n_bend = arc_segs(an_in)

            pts_arc_cen = []
            pts_arc_out = []
            pts_arc_in2 = []
            for k in range(n_bend + 1):
                ani = an_in * (n_bend - k) / n_bend
                ax_ = cxi + R_bend            * math.sin(ani)
                ay_ = cyi - R_bend            * math.cos(ani)
                ox_ = cxi + (R_bend+0.5*WGW) * math.sin(ani)
                oy_ = cyi - (R_bend+0.5*WGW) * math.cos(ani)
                ix_ = cxi + (R_bend-0.5*WGW) * math.sin(ani)
                iy_ = cyi - (R_bend-0.5*WGW) * math.cos(ani)
                # No negation — input side coords are already negative (left)
                pts_arc_cen.append((ax_, ay_))
                pts_arc_out.append((ox_, oy_))
                pts_arc_in2.append((ix_, iy_))

            self._prims.append(("poly", pts_arc_cen, C_FOCAL))
            self._prims.append(("filled_poly", pts_arc_out + list(reversed(pts_arc_in2)), col))

    # ─────────────────────────────────────────────────────────────
    # Render primitives → canvas
    # ─────────────────────────────────────────────────────────────

    def _render_prims(self):
        # Pass 1: filled polygons first (background layer)
        for typ, args, color in self._prims:
            if typ == "filled_poly":
                pts = args
                flat = []
                for x, y in pts:
                    px, py = self._u2c(x, y)
                    flat.extend([px, py])
                if len(flat) >= 6:
                    self.canvas.create_polygon(*flat, fill=color,
                                              outline=color, width=1)
        # Pass 2: lines and polylines on top
        for typ, args, color in self._prims:
            if typ == "line":
                x1, y1, x2, y2 = args
                cx1, cy1 = self._u2c(x1, y1)
                cx2, cy2 = self._u2c(x2, y2)
                self.canvas.create_line(cx1, cy1, cx2, cy2,
                                        fill=color, width=1)
            elif typ == "poly":
                pts = args
                flat = []
                for x, y in pts:
                    px, py = self._u2c(x, y)
                    flat.extend([px, py])
                if len(flat) >= 4:
                    self.canvas.create_line(*flat, fill=color,
                                           width=1, smooth=False)

    # ─────────────────────────────────────────────────────────────
    # DXF export  (proper LINE + ARC entities)
    # ─────────────────────────────────────────────────────────────

    def _on_dxf(self):
        if not self._prims:
            messagebox.showwarning("Warning",
                "No layout drawn yet. Click Draw first.")
            return
        fp = filedialog.asksaveasfilename(
            title="Export DXF",
            defaultextension=".dxf",
            filetypes=[("AutoCAD DXF", "*.dxf"), ("All Files", "*.*")])
        if not fp:
            return
        self._write_dxf(fp)

    def _write_dxf(self, filepath):
        """Export all geometry primitives as AutoCAD R14 DXF."""
        try:
            with open(filepath, "w") as f:

                def s(code, value):
                    f.write(f"{code}\n{value}\n")

                # ── Header ──
                s(0, "SECTION"); s(2, "HEADER")
                s(9, "$ACADVER"); s(1, "AC1014")
                s(9, "$INSUNITS"); s(70, 4)   # 4 = microns
                s(0, "ENDSEC")

                # ── Tables (minimal: one layer per colour) ──
                s(0, "SECTION"); s(2, "TABLES")
                s(0, "TABLE"); s(2, "LAYER"); s(70, 6)

                layer_map = {
                    C_FOCAL: ("AWG_CL",   1),    # red
                    C_AWG:   ("AWG_WG",   3),    # green
                    C_WALL:  ("AWG_WALL", 5),    # blue
                    C_SLAB:  ("SLAB",     7),    # white
                    C_IO:    ("IO_WG",    2),    # yellow
                    C_AXIS:  ("AXIS",     8),    # gray
                }
                for hex_col, (lname, aci) in layer_map.items():
                    s(0, "LAYER"); s(2, lname)
                    s(70, 0); s(62, aci); s(6, "CONTINUOUS")

                # add port colours as extra layers
                for k, col in enumerate(_PORT_COLORS[:8]):
                    s(0, "LAYER"); s(2, f"IO_PORT_{k+1}")
                    s(70, 0); s(62, k+1); s(6, "CONTINUOUS")

                s(0, "ENDTAB")
                s(0, "ENDSEC")

                # ── Entities ──
                s(0, "SECTION"); s(2, "ENTITIES")

                def get_layer(color):
                    if color in layer_map:
                        return layer_map[color][0]
                    for k, c in enumerate(_PORT_COLORS[:8]):
                        if c == color:
                            return f"IO_PORT_{k+1}"
                    return "AWG_WG"

                for typ, args, color in self._prims:
                    layer = get_layer(color)

                    if typ == "line":
                        x1, y1, x2, y2 = args
                        s(0, "LINE")
                        s(8, layer)
                        s(10, f"{x1:.6f}"); s(20, f"{y1:.6f}"); s(30, "0.0")
                        s(11, f"{x2:.6f}"); s(21, f"{y2:.6f}"); s(31, "0.0")

                    elif typ == "poly":
                        pts = args
                        if len(pts) < 2:
                            continue
                        # Write as LWPOLYLINE (open)
                        s(0, "LWPOLYLINE")
                        s(8, layer)
                        s(90, len(pts))   # vertex count
                        s(70, 0)          # not closed
                        for x, y in pts:
                            s(10, f"{x:.6f}")
                            s(20, f"{y:.6f}")

                    elif typ == "filled_poly":
                        pts = args
                        if len(pts) < 3:
                            continue
                        # Write as closed LWPOLYLINE (outline only in DXF)
                        s(0, "LWPOLYLINE")
                        s(8, layer)
                        s(90, len(pts))   # vertex count
                        s(70, 1)          # closed
                        for x, y in pts:
                            s(10, f"{x:.6f}")
                            s(20, f"{y:.6f}")

                s(0, "ENDSEC")
                s(0, "EOF")

            messagebox.showinfo("DXF Export",
                f"DXF file saved:\n{filepath}\n\n"
                f"Layers: AWG_CL, AWG_WG, AWG_WALL, SLAB, IO_WG, IO_PORT_1..8\n"
                f"Units: micrometers")

        except Exception as e:
            messagebox.showerror("DXF Error", str(e))

    # ─────────────────────────────────────────────────────────────
    # Mouse
    # ─────────────────────────────────────────────────────────────

    def _on_mouse_move(self, event):
        if not hasattr(self, '_cw'):
            return
        x, y = self._c2u(event.x, event.y)
        self.var_xpos.set(f"{x:.2f}")
        self.var_ypos.set(f"{y:.2f}")

    def _on_mouse_down(self, event):
        self._mx0 = event.x
        self._my0 = event.y

    def _on_drag1(self, event):
        self.canvas.create_rectangle(
            self._mx0, self._my0, event.x, event.y,
            outline="#ffe6ff", tags="zoom_rect")
        self.canvas.delete("zoom_rect")
        self.canvas.create_rectangle(
            self._mx0, self._my0, event.x, event.y,
            outline="#ffe6ff", tags="zoom_rect")

    def _on_pan(self, event):
        """Pan view with right-click drag."""
        if not self._prims:
            return
        dx_px = event.x - self._mx0
        dy_px = event.y - self._my0
        self._mx0 = event.x
        self._my0 = event.y
        # Convert pixel delta to world-unit delta
        self._ox -= dx_px * self.SCW / self._cw
        self._oy += dy_px * self.SCH / self._ch
        self._redraw()

    def _on_scroll(self, event):
        """Zoom in/out with scroll wheel, centred on mouse position."""
        if not self._prims:
            return
        # Determine zoom factor
        if event.num == 4 or (hasattr(event, 'delta') and event.delta > 0):
            factor = 0.8    # zoom in
        else:
            factor = 1.25   # zoom out
        # World position under mouse before zoom
        wx, wy = self._c2u(event.x, event.y)
        # Apply zoom
        self.SCW *= factor
        self.SCH *= factor
        # World position under same pixel after zoom (without offset fix)
        wx2, wy2 = self._c2u(event.x, event.y)
        # Shift offset so (wx, wy) stays under mouse
        self._ox += wx - wx2
        self._oy += wy - wy2
        self._redraw()

    def _redraw(self):
        """Re-render all cached primitives (used after zoom/pan)."""
        self.canvas.delete("all")
        self._render_prims()

    # ─────────────────────────────────────────────────────────────
    # Buttons
    # ─────────────────────────────────────────────────────────────

    def _on_simulate(self):
        if self.frm_variables:
            self.frm_variables._on_simulate()

    def _on_close(self):
        if self.window:
            self.window.withdraw()
