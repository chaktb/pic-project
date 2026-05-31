"""
Layout Parameters Form
Ported from Frm_Plot.frm (VB 5.0/6.0) to Python/Tkinter
"""

import tkinter as tk
from tkinter import ttk, messagebox
import globals as g
import awg_calc


class FrmPlot:
    """Device layout parameter input window."""

    def __init__(self, master, frm_variables=None):
        self.master = master
        self.frm_variables = frm_variables
        self.window = None
        self._frm_layout = None

        # StringVars for file I/O
        self.var_total_l   = tk.StringVar()
        self.var_d_io_port = tk.StringVar()
        self.var_d_slabs   = tk.StringVar()
        self.var_slb_angle = tk.StringVar()
        self.var_l_min     = tk.StringVar()
        self.var_b_min     = tk.StringVar()
        self.var_dum_i     = tk.StringVar(value="0")
        self.var_dum_o     = tk.StringVar(value="0")
        self.var_porder    = tk.StringVar()
        self.var_gap_sw    = tk.StringVar(value="0")
        self.var_fsr       = tk.StringVar()
        self.var_rowland_I = tk.IntVar(value=1)
        self.var_rowland_O = tk.IntVar(value=1)

        # Read-only calculated vars
        self.var_b_min_c   = tk.StringVar()
        self.var_l_min_c   = tk.StringVar()
        self.var_a_avg     = tk.StringVar()

        # Bind traces
        self.var_porder.trace_add("write", lambda *a: self._on_porder_change())

    def show(self):
        if self.window is None or not self.window.winfo_exists():
            self._build()
        self.window.deiconify()
        self.window.lift()

    def _build(self):
        self.window = tk.Toplevel(self.master)
        self.window.title("Plot Parameters")
        self.window.resizable(False, False)
        self.window.protocol("WM_DELETE_WINDOW", self._on_close)
        self._build_ui()

    def _build_ui(self):
        frm = ttk.LabelFrame(self.window,
                              text="Basic Structure Parameters", padding=6)
        frm.grid(row=0, column=0, padx=8, pady=6, sticky="nsew")

        def row(label, var, unit, readonly=False, r=None):
            idx = row.counter
            row.counter += 1
            ttk.Label(frm, text=label, anchor="w")\
                .grid(row=idx, column=0, sticky="w", padx=4, pady=2)
            state = "readonly" if readonly else "normal"
            e = ttk.Entry(frm, textvariable=var, width=10,
                          state=state, justify="right")
            e.grid(row=idx, column=1, padx=4, pady=2, sticky="w")
            if unit:
                ttk.Label(frm, text=unit).grid(row=idx, column=2, sticky="w")
            return e

        row.counter = 0

        row("total length of device",         self.var_total_l,   "um")
        row("io port distance",                self.var_d_io_port, "um")
        row("distance between slabs",          self.var_d_slabs,   "um")
        row("slab angle",                      self.var_slb_angle, "degree")
        row("minimum length of straight line", self.var_l_min,     "um")
        row("minium bend radius",              self.var_b_min,     "um")

        # Dummy waveguides (2 side-by-side entries)
        idx = row.counter
        row.counter += 1
        ttk.Label(frm, text="number of dummy waveguide")\
            .grid(row=idx, column=0, sticky="w", padx=4, pady=2)
        ttk.Entry(frm, textvariable=self.var_dum_i, width=4, justify="right")\
            .grid(row=idx, column=1, sticky="w", padx=4)
        ttk.Entry(frm, textvariable=self.var_dum_o, width=4, justify="right")\
            .grid(row=idx, column=2, sticky="w")

        # Phase order + FSR
        idx = row.counter
        row.counter += 1
        ttk.Label(frm, text="Phase Order")\
            .grid(row=idx, column=0, sticky="w", padx=4, pady=2)
        ttk.Entry(frm, textvariable=self.var_porder, width=7, justify="right")\
            .grid(row=idx, column=1, sticky="w", padx=4)
        ttk.Label(frm, text="FSR:")\
            .grid(row=idx, column=2, sticky="w")
        ttk.Entry(frm, textvariable=self.var_fsr, width=9, state="readonly",
                  justify="right")\
            .grid(row=idx, column=3, sticky="w")
        ttk.Label(frm, text="nm").grid(row=idx, column=4, sticky="w")

        row("gap between slab – waveguide", self.var_gap_sw, "um")

        row("awg length average",      self.var_a_avg,   "um", readonly=True)

        row("minium bend radius calc.", self.var_b_min_c, "um", readonly=True)
        row("minium length calc.",      self.var_l_min_c, "um", readonly=True)

        # Rowland circle checkboxes
        rl_frm = ttk.LabelFrame(self.window, text="Rowland Circle", padding=4)
        rl_frm.grid(row=1, column=0, padx=8, pady=4, sticky="ew")
        ttk.Checkbutton(rl_frm, text="Input Side",
                        variable=self.var_rowland_I)\
            .grid(row=0, column=0, padx=8)
        ttk.Checkbutton(rl_frm, text="Output Side",
                        variable=self.var_rowland_O)\
            .grid(row=0, column=1, padx=8)

        # Bottom buttons
        btn_frm = ttk.Frame(self.window)
        btn_frm.grid(row=2, column=0, pady=6, padx=8, sticky="ew")

        ttk.Button(btn_frm, text="Add Items", command=self._on_add_items)\
            .pack(side=tk.LEFT, padx=4)
        ttk.Button(btn_frm, text="Draw",      command=self._on_draw)\
            .pack(side=tk.LEFT, padx=4)
        ttk.Button(btn_frm, text="D Adjust",  command=self._on_d_adjust)\
            .pack(side=tk.LEFT, padx=4)
        self.btn_stop = ttk.Button(btn_frm, text="Stop",
                                   command=self._on_stop)
        self.btn_stop.pack(side=tk.LEFT, padx=4)
        ttk.Button(btn_frm, text="Close",     command=self._on_close)\
            .pack(side=tk.RIGHT, padx=4)

        self._is_adjusting = False

    # ─────────────────────────────────────────────────────────────
    # Handlers
    # ─────────────────────────────────────────────────────────────

    def _on_porder_change(self):
        try:
            mphase = float(self.var_porder.get())
            if mphase <= 0:
                return
            g.mphase = mphase
            if g.n_core * g.n_over * g.n_buff == 0:
                return
            awg_calc.m_phase(mphase)
            self.var_fsr.set(f"{g.fsr * 1e9:.4f}")
        except (ValueError, ZeroDivisionError):
            pass

    def _on_draw(self):
        frm_layout = self._get_frm_layout()
        frm_layout.draw()

    def _on_add_items(self):
        from frm_add import FrmAdd
        if not hasattr(self, '_frm_add') or self._frm_add is None:
            self._frm_add = FrmAdd(self.master)
        self._frm_add.show()

    def _on_d_adjust(self):
        """Auto-adjust slab distance to meet bend radius constraint."""
        self._is_adjusting = True
        delta = 1000.0

        self._do_draw_and_adjust(delta)

    def _do_draw_and_adjust(self, delta):
        if not self._is_adjusting:
            return
        frm_layout = self._get_frm_layout()
        frm_layout.draw()

        try:
            b_min_c = float(self.var_b_min_c.get())
            b_min   = float(self.var_b_min.get())
        except ValueError:
            self._is_adjusting = False
            return

        if abs(b_min_c - b_min) < 0.1 or delta < 1.0:
            self._is_adjusting = False
            return

        try:
            d = float(self.var_d_slabs.get())
        except ValueError:
            self._is_adjusting = False
            return

        if b_min_c > b_min:
            self.var_d_slabs.set(f"{d - delta:.2f}")
        else:
            self.var_d_slabs.set(f"{d + delta:.2f}")

        # Re-check after small delay
        self.window.after(100, lambda: self._do_draw_and_adjust(delta / 3.0))

    def _on_stop(self):
        self._is_adjusting = False

    def _on_close(self):
        self._is_adjusting = False
        if self.window:
            self.window.withdraw()

    def _get_frm_layout(self):
        if self._frm_layout is None:
            from frm_layout import FrmLayout
            self._frm_layout = FrmLayout(self.master, self.frm_variables, self)
        return self._frm_layout

    # ─────────────────────────────────────────────────────────────
    # Helpers used by layout form to update calculated fields
    # ─────────────────────────────────────────────────────────────

    def set_b_min_c(self, value, warn=False):
        self.var_b_min_c.set(f"{value:.2f}")

    def set_l_min_c(self, value, warn=False):
        self.var_l_min_c.set(f"{value:.5f}")

    def set_a_avg(self, value):
        self.var_a_avg.set(f"{value:.1f}")
