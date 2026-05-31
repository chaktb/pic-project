"""
Add Design Items Form (Tapers / Parabolic Waveguides / Segments)
Ported from Frm_Add.frm (VB 5.0/6.0) to Python/Tkinter

Values are written to globals immediately on every change.
"""

import tkinter as tk
from tkinter import ttk
import math
import globals as g


class FrmAdd:
    """Add tapers, parabolic waveguides, and segments to inputs/outputs."""

    def __init__(self, master):
        self.master = master
        self.window = None

        # Alpha display — created here so traces don't crash before _build_ui
        self.var_alpha = tk.StringVar(value="")

        # ── Input Port ──
        self.var_in_taper  = tk.IntVar(value=g.nchk_in_t)
        self.var_in_parab  = tk.IntVar(value=g.nchk_in_p)
        self.var_in_w      = tk.StringVar(value=str(g.W_in_item  * 1e6))
        self.var_in_l      = tk.StringVar(value=str(g.L_in_item  * 1e6))
        self.var_in_seg    = tk.IntVar(value=g.nchk_in_s)
        self.var_in_nseg   = tk.StringVar(value=str(g.N_in_seg))
        self.var_in_por    = tk.StringVar(value=str(g.Por_in_seg))

        # ── Output Port ──
        self.var_out_taper = tk.IntVar(value=g.nchk_out_t)
        self.var_out_parab = tk.IntVar(value=g.nchk_out_p)
        self.var_out_w     = tk.StringVar(value=str(g.W_out_item * 1e6))
        self.var_out_l     = tk.StringVar(value=str(g.L_out_item * 1e6))
        self.var_out_seg   = tk.IntVar(value=g.nchk_out_s)
        self.var_out_nseg  = tk.StringVar(value=str(g.N_out_seg))
        self.var_out_por   = tk.StringVar(value=str(g.Por_out_seg))

        # ── Input AWG side ──
        self.var_sin_taper = tk.IntVar(value=g.nchk_Sin_t)
        self.var_sin_parab = tk.IntVar(value=g.nchk_Sin_p)
        self.var_sin_w     = tk.StringVar(value=str(g.W_Sin_item * 1e6))
        self.var_sin_l     = tk.StringVar(value=str(g.L_Sin_item * 1e6))
        self.var_sin_seg   = tk.IntVar(value=g.nchk_Sin_s)
        self.var_sin_nseg  = tk.StringVar(value=str(g.N_Sin_seg))
        self.var_sin_por   = tk.StringVar(value=str(g.Por_Sin_seg))

        # ── Output AWG side ──
        self.var_sout_taper = tk.IntVar(value=g.nchk_Sout_t)
        self.var_sout_parab = tk.IntVar(value=g.nchk_Sout_p)
        self.var_sout_w     = tk.StringVar(value=str(g.W_Sout_item * 1e6))
        self.var_sout_l     = tk.StringVar(value=str(g.L_Sout_item * 1e6))
        self.var_sout_seg   = tk.IntVar(value=g.nchk_Sout_s)
        self.var_sout_nseg  = tk.StringVar(value=str(g.N_Sout_seg))
        self.var_sout_por   = tk.StringVar(value=str(g.Por_Sout_seg))

        # Register traces — every change writes to globals immediately
        for var in (self.var_in_w,   self.var_in_l,   self.var_in_nseg,   self.var_in_por,
                    self.var_out_w,  self.var_out_l,  self.var_out_nseg,  self.var_out_por,
                    self.var_sin_w,  self.var_sin_l,  self.var_sin_nseg,  self.var_sin_por,
                    self.var_sout_w, self.var_sout_l, self.var_sout_nseg, self.var_sout_por):
            var.trace_add("write", self._on_change)
        for var in (self.var_in_taper,   self.var_in_parab,   self.var_in_seg,
                    self.var_out_taper,  self.var_out_parab,  self.var_out_seg,
                    self.var_sin_taper,  self.var_sin_parab,  self.var_sin_seg,
                    self.var_sout_taper, self.var_sout_parab, self.var_sout_seg):
            var.trace_add("write", self._on_change)

    def show(self):
        if self.window is None or not self.window.winfo_exists():
            self._build()
        self._load_from_globals()
        self.window.deiconify()
        self.window.lift()

    def _build(self):
        self.window = tk.Toplevel(self.master)
        self.window.title("Add Items – Tapers / Parabolic / Segments")
        self.window.resizable(False, False)
        self.window.protocol("WM_DELETE_WINDOW", self._on_close)
        self._build_ui()

    def _build_ui(self):
        sections = [
            ("Input Port",
             self.var_in_taper,  self.var_in_parab,
             self.var_in_w,      self.var_in_l,
             self.var_in_seg,    self.var_in_nseg,  self.var_in_por),
            ("Output Port",
             self.var_out_taper, self.var_out_parab,
             self.var_out_w,     self.var_out_l,
             self.var_out_seg,   self.var_out_nseg, self.var_out_por),
            ("Input AWG side",
             self.var_sin_taper, self.var_sin_parab,
             self.var_sin_w,     self.var_sin_l,
             self.var_sin_seg,   self.var_sin_nseg, self.var_sin_por),
            ("Output AWG side",
             self.var_sout_taper, self.var_sout_parab,
             self.var_sout_w,     self.var_sout_l,
             self.var_sout_seg,   self.var_sout_nseg, self.var_sout_por),
        ]

        for col, (title, vt, vp, vw, vl, vs, vns, vpor) in enumerate(sections):
            frm = ttk.LabelFrame(self.window, text=title, padding=6)
            frm.grid(row=0, column=col, padx=6, pady=4, sticky="nsew")

            # ① Taper / Parabolic checkboxes
            ttk.Checkbutton(frm, text="Tapered Waveguide",  variable=vt)\
                .grid(row=0, column=0, columnspan=2, sticky="w")
            ttk.Checkbutton(frm, text="Parabolic Waveguide", variable=vp)\
                .grid(row=1, column=0, columnspan=2, sticky="w")

            # ② Width
            ttk.Label(frm, text="Width").grid(row=2, column=0, sticky="w")
            ttk.Entry(frm, textvariable=vw, width=8).grid(row=2, column=1, sticky="w")
            ttk.Label(frm, text="um").grid(row=2, column=2, sticky="w")

            # ③ Length
            ttk.Label(frm, text="Length").grid(row=3, column=0, sticky="w")
            ttk.Entry(frm, textvariable=vl, width=8).grid(row=3, column=1, sticky="w")
            ttk.Label(frm, text="um").grid(row=3, column=2, sticky="w")

            # ④ Alpha (read-only)
            ttk.Label(frm, text="alpha").grid(row=4, column=0, sticky="w")
            ttk.Entry(frm, textvariable=self.var_alpha, width=8, state="readonly")\
                .grid(row=4, column=1, sticky="w")

            ttk.Separator(frm, orient="horizontal")\
                .grid(row=5, column=0, columnspan=3, sticky="ew", pady=4)

            # ⑤ Segment checkbox
            ttk.Checkbutton(frm, text="Segment", variable=vs)\
                .grid(row=6, column=0, columnspan=2, sticky="w")

            # ⑥ N segment
            ttk.Label(frm, text="N segment").grid(row=7, column=0, sticky="w")
            ttk.Entry(frm, textvariable=vns, width=6).grid(row=7, column=1, sticky="w")

            # ⑦ Portion (%)
            ttk.Label(frm, text="Portion").grid(row=8, column=0, sticky="w")
            ttk.Entry(frm, textvariable=vpor, width=6).grid(row=8, column=1, sticky="w")
            ttk.Label(frm, text="%").grid(row=8, column=2, sticky="w")

        ttk.Button(self.window, text="Close", command=self._on_close)\
            .grid(row=1, column=0, columnspan=4, pady=6)

        self._update_alpha()

    # ─────────────────────────────────────────────────────────────
    # Sync globals ↔ UI
    # ─────────────────────────────────────────────────────────────

    def _on_change(self, *_):
        self._apply()
        self._update_alpha()

    def _apply(self):
        """Write all form values to globals immediately."""
        g.nchk_in_t   = self.var_in_taper.get()
        g.nchk_in_p   = self.var_in_parab.get()
        g.W_in_item   = self._flt(self.var_in_w)  * 1e-6
        g.L_in_item   = self._flt(self.var_in_l)  * 1e-6
        g.nchk_in_s   = self.var_in_seg.get()
        g.N_in_seg    = max(1, int(self._flt(self.var_in_nseg)))
        g.Por_in_seg  = max(0.0, min(100.0, self._flt(self.var_in_por)))

        g.nchk_out_t  = self.var_out_taper.get()
        g.nchk_out_p  = self.var_out_parab.get()
        g.W_out_item  = self._flt(self.var_out_w) * 1e-6
        g.L_out_item  = self._flt(self.var_out_l) * 1e-6
        g.nchk_out_s  = self.var_out_seg.get()
        g.N_out_seg   = max(1, int(self._flt(self.var_out_nseg)))
        g.Por_out_seg = max(0.0, min(100.0, self._flt(self.var_out_por)))

        g.nchk_Sin_t  = self.var_sin_taper.get()
        g.nchk_Sin_p  = self.var_sin_parab.get()
        g.W_Sin_item  = self._flt(self.var_sin_w) * 1e-6
        g.L_Sin_item  = self._flt(self.var_sin_l) * 1e-6
        g.nchk_Sin_s  = self.var_sin_seg.get()
        g.N_Sin_seg   = max(1, int(self._flt(self.var_sin_nseg)))
        g.Por_Sin_seg = max(0.0, min(100.0, self._flt(self.var_sin_por)))

        g.nchk_Sout_t  = self.var_sout_taper.get()
        g.nchk_Sout_p  = self.var_sout_parab.get()
        g.W_Sout_item  = self._flt(self.var_sout_w) * 1e-6
        g.L_Sout_item  = self._flt(self.var_sout_l) * 1e-6
        g.nchk_Sout_s  = self.var_sout_seg.get()
        g.N_Sout_seg   = max(1, int(self._flt(self.var_sout_nseg)))
        g.Por_Sout_seg = max(0.0, min(100.0, self._flt(self.var_sout_por)))

        g.f_saved = g.NO_

    def _load_from_globals(self):
        self.var_in_taper.set(g.nchk_in_t)
        self.var_in_parab.set(g.nchk_in_p)
        self.var_in_w.set(str(g.W_in_item  * 1e6))
        self.var_in_l.set(str(g.L_in_item  * 1e6))
        self.var_in_seg.set(g.nchk_in_s)
        self.var_in_nseg.set(str(g.N_in_seg))
        self.var_in_por.set(str(g.Por_in_seg))

        self.var_out_taper.set(g.nchk_out_t)
        self.var_out_parab.set(g.nchk_out_p)
        self.var_out_w.set(str(g.W_out_item * 1e6))
        self.var_out_l.set(str(g.L_out_item * 1e6))
        self.var_out_seg.set(g.nchk_out_s)
        self.var_out_nseg.set(str(g.N_out_seg))
        self.var_out_por.set(str(g.Por_out_seg))

        self.var_sin_taper.set(g.nchk_Sin_t)
        self.var_sin_parab.set(g.nchk_Sin_p)
        self.var_sin_w.set(str(g.W_Sin_item * 1e6))
        self.var_sin_l.set(str(g.L_Sin_item * 1e6))
        self.var_sin_seg.set(g.nchk_Sin_s)
        self.var_sin_nseg.set(str(g.N_Sin_seg))
        self.var_sin_por.set(str(g.Por_Sin_seg))

        self.var_sout_taper.set(g.nchk_Sout_t)
        self.var_sout_parab.set(g.nchk_Sout_p)
        self.var_sout_w.set(str(g.W_Sout_item * 1e6))
        self.var_sout_l.set(str(g.L_Sout_item * 1e6))
        self.var_sout_seg.set(g.nchk_Sout_s)
        self.var_sout_nseg.set(str(g.N_Sout_seg))
        self.var_sout_por.set(str(g.Por_Sout_seg))

    def _flt(self, var):
        try:
            return float(var.get())
        except (ValueError, AttributeError):
            return 0.0

    def _update_alpha(self):
        try:
            w_um  = self._flt(self.var_in_w)
            l_um  = self._flt(self.var_in_l)
            wgw_um = g.wg_w * 1e6
            if l_um > 0 and w_um > wgw_um:
                alpha = math.degrees(math.atan((w_um - wgw_um) / (2.0 * l_um)))
                self.var_alpha.set(f"{alpha:.4f} deg")
            else:
                self.var_alpha.set("")
        except (ValueError, ZeroDivisionError):
            self.var_alpha.set("")

    def _on_close(self):
        self._apply()
        if self.window:
            self.window.withdraw()
