"""
AWG WDM Design Utility - Main Design Entry Form
Ported from Frm_Variables.frm (VB 5.0/6.0) to Python/Tkinter

Corresponds to: Frm_Variables (startup form)
"""

import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import globals as g
import awg_calc


class FrmVariables:
    """Main AWG design parameter input window."""

    def __init__(self, master):
        self.master = master

        # Use the root Tk window directly (ensures menu bar works on all platforms)
        self.window = master
        self.window.title(g.VER_STRING)
        self.window.resizable(False, False)
        self.window.protocol("WM_DELETE_WINDOW", self._on_exit)

        # References to other forms (created lazily)
        self._frm_param    = None
        self._frm_plot     = None
        self._frm_layout   = None
        self._frm_index    = None
        self._frm_about    = None
        self._frm_add_items = None

        self._build_menu()
        self._build_ui()
        self._reset_all()

    # ─────────────────────────────────────────────────────────────
    # Menu
    # ─────────────────────────────────────────────────────────────

    def _build_menu(self):
        menubar = tk.Menu(self.window)

        file_menu = tk.Menu(menubar, tearoff=0)
        file_menu.add_command(label="New",              command=self._on_new)
        file_menu.add_command(label="Open",             command=self._on_open)
        file_menu.add_command(label="Open Old Version", command=self._on_open_old)
        file_menu.add_command(label="Save As",          command=self._on_save_as)
        file_menu.add_separator()
        file_menu.add_command(label="Exit",             command=self._on_exit)
        menubar.add_cascade(label="File", menu=file_menu)

        wg_menu = tk.Menu(menubar, tearoff=0)
        wg_menu.add_command(label="Refractive Index Calc...", command=self._on_waveguide)
        menubar.add_cascade(label="Waveguide", menu=wg_menu)

        help_menu = tk.Menu(menubar, tearoff=0)
        help_menu.add_command(label="About", command=self._on_about)
        menubar.add_cascade(label="Help", menu=help_menu)

        self.window.config(menu=menubar)

    # ─────────────────────────────────────────────────────────────
    # UI Layout
    # ─────────────────────────────────────────────────────────────

    def _build_ui(self):
        pad = {"padx": 5, "pady": 3}

        # ── Wavelength Frame ──
        frm_wl = ttk.LabelFrame(self.window, text="Wavelength")
        frm_wl.grid(row=0, column=0, padx=8, pady=6, sticky="nsew")
        self._build_wavelength_frame(frm_wl)

        # ── Waveguide Basics Frame ──
        frm_wg = ttk.LabelFrame(self.window, text="Waveguide Basics")
        frm_wg.grid(row=0, column=1, padx=8, pady=6, sticky="nsew")
        self._build_waveguide_frame(frm_wg)

        # ── AWG Frame ──
        frm_awg = ttk.LabelFrame(self.window, text="Array Waveguide Gratings")
        frm_awg.grid(row=0, column=2, padx=8, pady=6, sticky="nsew")
        self._build_awg_frame(frm_awg)

        # ── Bottom Buttons ──
        btn_frame = ttk.Frame(self.window)
        btn_frame.grid(row=1, column=0, columnspan=3, pady=6, sticky="e", padx=8)

        ttk.Button(btn_frame, text="Simulate", width=10, command=self._on_simulate)\
            .pack(side=tk.LEFT, padx=4)
        ttk.Button(btn_frame, text="Draw",     width=10, command=self._on_draw)\
            .pack(side=tk.LEFT, padx=4)
        ttk.Button(btn_frame, text="Exit",     width=10, command=self._on_exit)\
            .pack(side=tk.LEFT, padx=4)

    # ─── Wavelength sub-frame ───────────────────────────────────

    def _build_wavelength_frame(self, parent):
        # Number of channels + 1xN / 2xN
        ttk.Label(parent, text="Number of Channels").grid(row=0, column=0, sticky="w", padx=4, pady=2)
        self.var_nc = tk.StringVar(value="2")
        nc_entry = ttk.Entry(parent, textvariable=self.var_nc, width=5)
        nc_entry.grid(row=0, column=1, sticky="w", pady=2)
        nc_entry.bind("<FocusOut>", lambda e: self._on_nc_change())

        self.var_1n = tk.IntVar(value=1)
        self.var_2n = tk.IntVar(value=0)
        ttk.Checkbutton(parent, text="1xN", variable=self.var_1n, command=self._on_1n_click)\
            .grid(row=0, column=2, padx=2)
        ttk.Checkbutton(parent, text="2xN", variable=self.var_2n, command=self._on_2n_click)\
            .grid(row=0, column=3, padx=2)

        # Center wavelength
        ttk.Label(parent, text="Center Wavelength").grid(row=1, column=0, sticky="w", padx=4, pady=2)
        self.var_wlc = tk.StringVar()
        wlc_entry = ttk.Entry(parent, textvariable=self.var_wlc, width=8)
        wlc_entry.grid(row=1, column=1, sticky="w", pady=2)
        wlc_entry.bind("<FocusOut>", lambda e: self._calc_wl())
        wlc_entry.bind("<Return>",   lambda e: self._calc_wl())
        ttk.Label(parent, text="nm").grid(row=1, column=2, sticky="w")

        # Channel spacing
        ttk.Label(parent, text="Channel Spacing").grid(row=2, column=0, sticky="w", padx=4, pady=2)
        self.var_wld = tk.StringVar()
        wld_entry = ttk.Entry(parent, textvariable=self.var_wld, width=8)
        wld_entry.grid(row=2, column=1, sticky="w", pady=2)
        wld_entry.bind("<FocusOut>", lambda e: self._calc_wl())
        wld_entry.bind("<Return>",   lambda e: self._calc_wl())
        ttk.Label(parent, text="nm").grid(row=2, column=2, sticky="w")

        # Start wavelength (read-only)
        ttk.Label(parent, text="Start Wavelength").grid(row=3, column=0, sticky="w", padx=4, pady=2)
        self.var_wls = tk.StringVar()
        ttk.Entry(parent, textvariable=self.var_wls, width=8, state="readonly")\
            .grid(row=3, column=1, sticky="w", pady=2)
        ttk.Label(parent, text="nm").grid(row=3, column=2, sticky="w")

        # End wavelength (read-only)
        ttk.Label(parent, text="End Wavelength").grid(row=4, column=0, sticky="w", padx=4, pady=2)
        self.var_wle = tk.StringVar()
        ttk.Entry(parent, textvariable=self.var_wle, width=8, state="readonly")\
            .grid(row=4, column=1, sticky="w", pady=2)
        ttk.Label(parent, text="nm").grid(row=4, column=2, sticky="w")

    # ─── Waveguide sub-frame ────────────────────────────────────

    def _build_waveguide_frame(self, parent):
        # All fields are editable.
        # Core / Over / Buff indices: also settable via Waveguide menu.
        rows = [
            ("Waveguide Width",        "var_w",    "um"),
            ("Waveguide Height",       "var_h",    "um"),
            ("Polarization (X/Y/S)",   "var_pol",  ""),
            ("Core Refractive Index",  "var_core", ""),
            ("Over Refractive Index",  "var_over", ""),
            ("Buff Refractive Index",  "var_buff", ""),
        ]
        for r, (label, attr, unit) in enumerate(rows):
            ttk.Label(parent, text=label).grid(row=r, column=0, sticky="w", padx=4, pady=2)
            var = tk.StringVar()
            setattr(self, attr, var)
            entry = ttk.Entry(parent, textvariable=var, width=10)
            entry.grid(row=r, column=1, sticky="w", pady=2)
            if unit:
                ttk.Label(parent, text=unit).grid(row=r, column=2, sticky="w")

            # Bind change handlers
            if attr == "var_w":
                var.trace_add("write", lambda *a: self._on_w_change())
            elif attr == "var_h":
                var.trace_add("write", lambda *a: self._on_h_change())
            elif attr == "var_pol":
                var.trace_add("write", lambda *a: self._on_pol_change())
                entry.bind("<Double-Button-1>", self._on_pol_dblclick)
            elif attr == "var_core":
                var.trace_add("write", lambda *a: self._on_core_change())
            elif attr == "var_over":
                var.trace_add("write", lambda *a: self._on_over_change())
            elif attr == "var_buff":
                var.trace_add("write", lambda *a: self._on_buff_change())

    # ─── AWG sub-frame ──────────────────────────────────────────

    def _build_awg_frame(self, parent):
        ttk.Label(parent, text="Number of Array Waveguides")\
            .grid(row=0, column=0, sticky="w", padx=4, pady=2)
        self.var_na = tk.StringVar()
        na_entry = ttk.Entry(parent, textvariable=self.var_na, width=8)
        na_entry.grid(row=0, column=1, sticky="w", pady=2)
        na_entry.bind("<FocusOut>", lambda e: self._on_na_change())
        na_entry.bind("<Return>",   lambda e: self._on_na_change())

        ttk.Label(parent, text="Dist. between Array WGs")\
            .grid(row=1, column=0, sticky="w", padx=4, pady=2)
        self.var_d_aw = tk.StringVar()
        daw_e = ttk.Entry(parent, textvariable=self.var_d_aw, width=8)
        daw_e.grid(row=1, column=1, sticky="w", pady=2)
        daw_e.bind("<FocusOut>", lambda e: self._on_daw_change())
        daw_e.bind("<Return>",   lambda e: self._on_daw_change())
        ttk.Label(parent, text="um").grid(row=1, column=2, sticky="w")

        ttk.Label(parent, text="Dist. between IO Waveguides")\
            .grid(row=2, column=0, sticky="w", padx=4, pady=2)
        self.var_d_io = tk.StringVar()
        dio_e = ttk.Entry(parent, textvariable=self.var_d_io, width=8)
        dio_e.grid(row=2, column=1, sticky="w", pady=2)
        dio_e.bind("<FocusOut>", lambda e: self._on_dio_change())
        dio_e.bind("<Return>",   lambda e: self._on_dio_change())
        ttk.Label(parent, text="um").grid(row=2, column=2, sticky="w")

    # ─────────────────────────────────────────────────────────────
    # Widget → Global variable handlers
    # ─────────────────────────────────────────────────────────────

    def _on_nc_change(self):
        try:
            nc = int(float(self.var_nc.get()))
            g.nc = max(nc, 1)
        except ValueError:
            return
        self._calc_wl()

    def _on_1n_click(self):
        g.f_saved = g.NO_
        g.f_completed = g.NO_
        if self.var_1n.get():
            self.var_2n.set(0)
        self._calc_wl()

    def _on_2n_click(self):
        g.f_saved = g.NO_
        g.f_completed = g.NO_
        if self.var_2n.get():
            self.var_1n.set(0)
        self._calc_wl()

    def _on_w_change(self):
        try:
            g.wg_w = float(self.var_w.get()) * 1e-6
        except ValueError:
            pass

    def _on_h_change(self):
        try:
            g.wg_h = float(self.var_h.get()) * 1e-6
        except ValueError:
            pass

    def _on_core_change(self):
        try:
            g.n_core = float(self.var_core.get())
        except ValueError:
            pass

    def _on_over_change(self):
        try:
            g.n_over = float(self.var_over.get())
        except ValueError:
            pass

    def _on_buff_change(self):
        try:
            g.n_buff = float(self.var_buff.get())
        except ValueError:
            pass

    def _on_pol_change(self):
        pol = self.var_pol.get().upper()
        if pol in ("", "X", "Y", "S"):
            g.pol = pol
        else:
            messagebox.showwarning("Warning", "Polarization must be X, Y, or S")
            self.var_pol.set("")

    def _on_pol_dblclick(self, event=None):
        pol = self.var_pol.get()
        if pol == "X":
            self.var_pol.set("Y")
        elif pol == "Y":
            self.var_pol.set("S")
        elif pol == "S":
            self.var_pol.set("X")
        else:
            self.var_pol.set("X")

    def _on_na_change(self):
        try:
            na = int(float(self.var_na.get()))
            g.na = max(na, 0)
        except ValueError:
            return
        g.f_saved = g.NO_
        g.f_completed = g.NO_
        g.init_arrays(g.nc, g.na)

    def _on_daw_change(self):
        try:
            g.d_aw = float(self.var_d_aw.get()) * 1e-6
        except ValueError:
            pass
        g.f_saved = g.NO_
        g.f_completed = g.NO_

    def _on_dio_change(self):
        try:
            g.d_io = float(self.var_d_io.get()) * 1e-6
        except ValueError:
            pass
        g.f_saved = g.NO_
        g.f_completed = g.NO_

    def _calc_wl(self):
        """Calculate start/end wavelengths from center + spacing + channel count."""
        try:
            nc  = int(float(self.var_nc.get()))
            wlc = float(self.var_wlc.get()) * 1e-9
            wld = float(self.var_wld.get()) * 1e-9
        except ValueError:
            return

        if nc <= 0:
            messagebox.showwarning("Warning", "Number of channels must be > 0")
            return

        g.nc  = nc
        g.wl_c = wlc
        g.wl_d = wld

        is_1n = self.var_1n.get()

        if (nc % 2 == 0) and not is_1n:
            wl_s = wlc - 0.5 * nc * wld
            wl_e = wlc + 0.5 * (nc - 2) * wld
        else:
            wl_s = wlc - 0.5 * (nc - 1) * wld
            wl_e = wlc + 0.5 * (nc - 1) * wld

        g.wl_s = wl_s
        g.wl_e = wl_e

        self.var_wls.set(f"{wl_s * 1e9:.4f}")
        self.var_wle.set(f"{wl_e * 1e9:.4f}")

    def _read_waveguide_params(self):
        """Sync all waveguide parameter widgets → globals."""
        try:
            g.nc  = int(float(self.var_nc.get()))
            g.wg_w = float(self.var_w.get()) * 1e-6
            g.wg_h = float(self.var_h.get()) * 1e-6
            g.pol  = self.var_pol.get().upper() or "X"
            g.n_core = float(self.var_core.get())
            g.n_over = float(self.var_over.get())
            g.n_buff = float(self.var_buff.get())
            g.na  = int(float(self.var_na.get()))
            g.d_aw = float(self.var_d_aw.get()) * 1e-6
            g.d_io = float(self.var_d_io.get()) * 1e-6
        except ValueError as e:
            messagebox.showerror("Input Error", str(e))
            return False
        return True

    # ─────────────────────────────────────────────────────────────
    # Button actions
    # ─────────────────────────────────────────────────────────────

    def _on_simulate(self):
        if not self._validate_for_simulate():
            return
        if not self._read_waveguide_params():
            return

        try:
            self._calc_wl()
            awg_calc.calc_fsr()
            # Phase order from Frm_Plot
            frm_plot = self._get_frm_plot()
            mphase_str = frm_plot.var_porder.get()
            try:
                mphase_val = float(mphase_str) if mphase_str else 0.0
            except ValueError:
                mphase_val = 0.0

            g.mphase = mphase_val
            awg_calc.m_phase(g.mphase if g.mphase > 0 else g.mphase)
            frm_plot.var_fsr.set(f"{g.fsr * 1e9:.4f}")
        except Exception as e:
            messagebox.showerror("Calculation Error", str(e))
            return

        # Show parameters form
        frm_param = self._get_frm_param()
        frm_param.show()
        frm_param.refresh_params()

    def _on_draw(self):
        if not self._validate_for_simulate():
            return
        if not self._read_waveguide_params():
            return

        try:
            self._calc_wl()
            awg_calc.calc_fsr()
            frm_plot = self._get_frm_plot()
            mphase_str = frm_plot.var_porder.get()
            try:
                mphase_val = float(mphase_str) if mphase_str else 0.0
            except ValueError:
                mphase_val = 0.0
            g.mphase = mphase_val
            awg_calc.m_phase(g.mphase if g.mphase > 0 else g.mphase)
            frm_plot.var_fsr.set(f"{g.fsr * 1e9:.4f}")
        except Exception as e:
            messagebox.showerror("Calculation Error", str(e))
            return

        frm_param = self._get_frm_param()
        frm_param.show()
        frm_param.refresh_params()

        frm_plot = self._get_frm_plot()
        frm_plot.show()

    def _on_exit(self):
        if g.f_completed == g.YES and g.f_saved == g.NO_:
            res = messagebox.askyesnocancel(
                "Exit", "File not saved. Save it now?")
            if res is None:   # Cancel
                return
            if res:           # Yes
                self._on_save_as()

        self.master.destroy()

    def _validate_for_simulate(self):
        try:
            n_core = float(self.var_core.get())
            if n_core < 1:
                messagebox.showwarning("Warning", "Waveguide not defined (core index < 1)")
                return False
        except ValueError:
            messagebox.showwarning("Warning", "Waveguide not defined")
            return False
        return True

    # ─────────────────────────────────────────────────────────────
    # Menu actions
    # ─────────────────────────────────────────────────────────────

    def _on_new(self):
        if g.f_completed == g.YES and g.f_saved == g.NO_:
            if not messagebox.askokcancel("New", "File not saved! Continue?"):
                return
        self._reset_all()

    def _on_open(self):
        filepath = filedialog.askopenfilename(
            title="Open AWG Design File",
            filetypes=[("AWG Design", "*.awg"), ("All Files", "*.*")])
        if not filepath:
            return
        self._load_file(filepath)

    def _on_open_old(self):
        """Open an older version AWG file (extra field skipped)."""
        filepath = filedialog.askopenfilename(
            title="Open Old AWG File",
            filetypes=[("AWG Design", "*.awg"), ("All Files", "*.*")])
        if not filepath:
            return
        self._load_file(filepath, old_version=True)

    def _on_save_as(self):
        if g.f_completed == g.NO_:
            messagebox.showwarning("Warning", "Design not completed")
            return
        filepath = filedialog.asksaveasfilename(
            title="Save AWG Design",
            defaultextension=".awg",
            filetypes=[("AWG Design", "*.awg"), ("All Files", "*.*")])
        if not filepath:
            return
        self._save_file(filepath)

    def _on_waveguide(self):
        """Open refractive index / waveguide mode calculator."""
        frm_index = self._get_frm_index()
        frm_index.set_inputs(
            wl_nm  = self.var_wlc.get(),
            w_um   = self.var_w.get(),
            h_um   = self.var_h.get(),
            pol    = self.var_pol.get(),
            n_core = self.var_core.get(),
            n_over = self.var_over.get(),
            n_buff = self.var_buff.get(),
        )
        frm_index.show()

    def _on_about(self):
        frm_about = self._get_frm_about()
        frm_about.show()

    # ─────────────────────────────────────────────────────────────
    # File I/O  (matches VB sequential file format)
    # ─────────────────────────────────────────────────────────────

    def _load_file(self, filepath, old_version=False):
        """Load .awg design file."""
        try:
            with open(filepath, "r") as fh:
                lines = [ln.strip() for ln in fh.readlines()]

            def read():
                return lines.pop(0) if lines else ""

            self._reset_all()

            self.var_nc.set(read())
            self.var_1n.set(int(read()))
            self.var_2n.set(int(read()))
            self.var_wlc.set(read())
            self.var_wld.set(read())
            self.var_w.set(read())
            self.var_h.set(read())
            if old_version:
                read()  # skip index delta field from old versions
            self.var_pol.set(read())
            self.var_na.set(read())
            self.var_d_aw.set(read())
            self.var_d_io.set(read())

            fp = self._get_frm_plot()
            fp.var_rowland_I.set(int(read()))
            fp.var_rowland_O.set(int(read()))
            fp.var_total_l.set(read())
            fp.var_d_io_port.set(read())
            fp.var_d_slabs.set(read())
            fp.var_slb_angle.set(read())
            fp.var_l_min.set(read())
            fp.var_b_min.set(read())
            fp.var_dum_i.set(read())
            fp.var_dum_o.set(read())
            fp.var_porder.set(read())
            fp.var_gap_sw.set(read())

            fa = self._get_frm_add()
            fa.var_in_taper.set(int(read()))
            fa.var_in_parab.set(int(read()))
            fa.var_in_w.set(read())
            fa.var_in_l.set(read())
            fa.var_out_taper.set(int(read()))
            fa.var_out_parab.set(int(read()))
            fa.var_out_w.set(read())
            fa.var_out_l.set(read())
            fa.var_sin_taper.set(int(read()))
            fa.var_sin_parab.set(int(read()))
            fa.var_sin_w.set(read())
            fa.var_sin_l.set(read())
            fa.var_sout_taper.set(int(read()))
            fa.var_sout_parab.set(int(read()))
            fa.var_sout_w.set(read())
            fa.var_sout_l.set(read())

            fl = self._get_frm_layout()
            fl.var_offset1.set(read())
            fl.var_offset2[0].set(read())
            fl.var_offset2[1].set(read())
            fl.var_offset2[2].set(read())
            fl.var_focal_awg.set(read())
            fl.var_focal_io.set(read())

            fp2 = self._get_frm_param()
            fp2.var_text.set(read())

            if lines:
                fl.var_wsl.set(read())

            self._calc_wl()
            self._read_waveguide_params()

            g.f_saved = g.YES
            g.f_completed = g.NO_
            self.window.title(f"{g.VER_STRING}  –  {filepath}")

        except Exception as e:
            messagebox.showerror("File Open Error", str(e))

    def _save_file(self, filepath):
        """Save .awg design file."""
        try:
            fp  = self._get_frm_plot()
            fa  = self._get_frm_add()
            fl  = self._get_frm_layout()
            fpm = self._get_frm_param()

            lines = [
                self.var_nc.get(),
                str(self.var_1n.get()),
                str(self.var_2n.get()),
                self.var_wlc.get(),
                self.var_wld.get(),
                self.var_w.get(),
                self.var_h.get(),
                self.var_pol.get(),
                self.var_na.get(),
                self.var_d_aw.get(),
                self.var_d_io.get(),
                str(fp.var_rowland_I.get()),
                str(fp.var_rowland_O.get()),
                fp.var_total_l.get(),
                fp.var_d_io_port.get(),
                fp.var_d_slabs.get(),
                fp.var_slb_angle.get(),
                fp.var_l_min.get(),
                fp.var_b_min.get(),
                fp.var_dum_i.get(),
                fp.var_dum_o.get(),
                fp.var_porder.get(),
                fp.var_gap_sw.get(),
                str(fa.var_in_taper.get()),
                str(fa.var_in_parab.get()),
                fa.var_in_w.get(),
                fa.var_in_l.get(),
                str(fa.var_out_taper.get()),
                str(fa.var_out_parab.get()),
                fa.var_out_w.get(),
                fa.var_out_l.get(),
                str(fa.var_sin_taper.get()),
                str(fa.var_sin_parab.get()),
                fa.var_sin_w.get(),
                fa.var_sin_l.get(),
                str(fa.var_sout_taper.get()),
                str(fa.var_sout_parab.get()),
                fa.var_sout_w.get(),
                fa.var_sout_l.get(),
                fl.var_offset1.get(),
                fl.var_offset2[0].get(),
                fl.var_offset2[1].get(),
                fl.var_offset2[2].get(),
                fl.var_focal_awg.get(),
                fl.var_focal_io.get(),
                fpm.var_text.get(),
                fl.var_wsl.get(),
            ]

            with open(filepath, "w") as fh:
                fh.write("\n".join(lines) + "\n")

            g.f_saved = g.YES
            self.window.title(f"{g.VER_STRING}  –  {filepath}")

        except Exception as e:
            messagebox.showerror("Save Error", str(e))

    # ─────────────────────────────────────────────────────────────
    # Reset
    # ─────────────────────────────────────────────────────────────

    def _reset_all(self):
        """Reset all fields to defaults."""
        self.var_nc.set("1")
        self.var_1n.set(1)
        self.var_2n.set(0)
        self.var_wlc.set("")
        self.var_wld.set("")
        self.var_wls.set("")
        self.var_wle.set("")
        self.var_w.set("")
        self.var_h.set("")
        self.var_pol.set("X")
        self.var_core.set("")
        self.var_over.set("")
        self.var_buff.set("")
        self.var_na.set("")
        self.var_d_aw.set("")
        self.var_d_io.set("")

        # Reset globals
        g.nc = 1
        g.na = 0
        g.wl_c = g.wl_d = g.wl_s = g.wl_e = 0.0
        g.wg_w = g.wg_h = 0.0
        g.pol = "X"
        g.n_core = g.n_over = g.n_buff = 0.0
        g.mphase = 0.0
        g.d_aw = g.d_io = 0.0
        g.f_saved = g.YES
        g.f_completed = g.NO_

        self.window.title(g.VER_STRING)

    def update_index_display(self, n_core, n_over, n_buff):
        """Called from frm_index after calculation to update read-only fields."""
        self.var_core.set(f"{n_core:.7f}")
        self.var_over.set(f"{n_over:.7f}")
        self.var_buff.set(f"{n_buff:.7f}")
        g.n_core = n_core
        g.n_over = n_over
        g.n_buff = n_buff

    # ─────────────────────────────────────────────────────────────
    # Lazy form getters
    # ─────────────────────────────────────────────────────────────

    def _get_frm_param(self):
        if self._frm_param is None:
            from frm_param import FrmParam
            self._frm_param = FrmParam(self.master, self)
        return self._frm_param

    def _get_frm_plot(self):
        if self._frm_plot is None:
            from frm_plot import FrmPlot
            self._frm_plot = FrmPlot(self.master, self)
        return self._frm_plot

    def _get_frm_layout(self):
        if self._frm_layout is None:
            from frm_layout import FrmLayout
            self._frm_layout = FrmLayout(self.master, self)
        return self._frm_layout

    def _get_frm_index(self):
        if self._frm_index is None:
            from frm_index import FrmIndex
            self._frm_index = FrmIndex(self.master, self)
        return self._frm_index

    def _get_frm_about(self):
        if self._frm_about is None:
            from frm_about import FrmAbout
            self._frm_about = FrmAbout(self.master)
        return self._frm_about

    def _get_frm_add(self):
        if self._frm_add_items is None:
            from frm_add import FrmAdd
            self._frm_add_items = FrmAdd(self.master)
        return self._frm_add_items
