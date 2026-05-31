"""
Design Parameters & Simulation Form
Ported from Frm_Param.frm (VB 5.0/6.0) to Python/Tkinter
"""

import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import threading
import math
import globals as g
import awg_calc


class FrmParam:
    """Design parameters display and simulation control."""

    def __init__(self, master, frm_variables=None):
        self.master = master
        self.frm_variables = frm_variables
        self.window = None
        self._sim_running = False
        self._sim_stop = False

        # StringVar exposed to file I/O
        self.var_text = tk.StringVar()

    def show(self):
        if self.window is None or not self.window.winfo_exists():
            self._build()
        self.window.deiconify()
        self.window.lift()

    def _build(self):
        self.window = tk.Toplevel(self.master)
        self.window.title("Design Parameters")
        self.window.resizable(False, False)
        self.window.protocol("WM_DELETE_WINDOW", self._on_close)

        self._build_menu()
        self._build_ui()

    def _build_menu(self):
        menubar = tk.Menu(self.window)
        file_m = tk.Menu(menubar, tearoff=0)
        file_m.add_command(label="Print", command=self._on_print)
        menubar.add_cascade(label="File", menu=file_m)
        self.window.config(menu=menubar)

    def _build_ui(self):
        # ── Left: Parameter display ──
        left = ttk.Frame(self.window, padding=4)
        left.grid(row=0, column=0, sticky="nsew")

        self.txt_params = tk.Text(left, width=55, height=22,
                                  font=("Courier", 9), state="disabled")
        self.txt_params.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        sb = ttk.Scrollbar(left, command=self.txt_params.yview)
        sb.pack(side=tk.RIGHT, fill=tk.Y)
        self.txt_params.configure(yscrollcommand=sb.set)

        # ── Right: Controls ──
        right = ttk.Frame(self.window, padding=4)
        right.grid(row=0, column=1, sticky="nsew")

        # Coupling coefficients
        cc_frm = ttk.LabelFrame(right, text="Coupling Coefficients", padding=4)
        cc_frm.grid(row=0, column=0, sticky="ew", pady=4)
        ttk.Button(cc_frm, text="Uniform", width=9, command=self._on_uniform)\
            .grid(row=0, column=0, padx=3)
        ttk.Button(cc_frm, text="Gauss",   width=9, command=self._on_gauss)\
            .grid(row=0, column=1, padx=3)
        ttk.Button(cc_frm, text="BPM",     width=9, command=self._on_bpm)\
            .grid(row=0, column=2, padx=3)

        # Gaussian width sub-panel
        self.gauss_frm = ttk.LabelFrame(right, text="Gaussian Width", padding=4)
        self.gauss_frm.grid(row=1, column=0, sticky="ew", pady=2)
        self.var_gauss_w = tk.StringVar(value="0.35")
        ttk.Entry(self.gauss_frm, textvariable=self.var_gauss_w, width=7)\
            .grid(row=0, column=0)
        ttk.Button(self.gauss_frm, text="Apply", command=self._apply_gauss)\
            .grid(row=0, column=1, padx=4)
        self.gauss_frm.grid_remove()  # hidden by default

        # Simulation variables
        sim_frm = ttk.LabelFrame(right, text="Simulation Variables", padding=4)
        sim_frm.grid(row=2, column=0, sticky="ew", pady=4)

        rows_sv = [
            ("waveguide loss",   "var_s_loss",    "dB/cm", "0"),
            ("waveguide width",  "var_s_width",   "μm",    "0"),
            ("waveguide height", "var_s_height",  "μm",    "0"),
            ("Temperature",      "var_s_tempr",   "°C",    "25"),
            ("WL start",         "var_s_wlstart", "nm",    "0"),
            ("WL end",           "var_s_wlend",   "nm",    "0"),
            ("Points",           "var_s_points",  "",      "500"),
        ]
        for r, (label, attr, unit, default) in enumerate(rows_sv):
            ttk.Label(sim_frm, text=label, anchor="w")\
                .grid(row=r, column=0, sticky="w", padx=2)
            var = tk.StringVar(value=default)
            setattr(self, attr, var)
            ttk.Entry(sim_frm, textvariable=var, width=9, justify="right")\
                .grid(row=r, column=1, padx=2, pady=1)
            if unit:
                ttk.Label(sim_frm, text=unit).grid(row=r, column=2, sticky="w")

        # Input port selector
        chan_frm = ttk.LabelFrame(right, text="Channel Selection", padding=4)
        chan_frm.grid(row=3, column=0, sticky="ew", pady=4)
        ttk.Label(chan_frm, text="Input port:").grid(row=0, column=0)
        self.var_input_port = tk.StringVar(value="1")
        ttk.Entry(chan_frm, textvariable=self.var_input_port, width=4)\
            .grid(row=0, column=1)
        ttk.Label(chan_frm, text="Output ports:").grid(row=1, column=0)
        self.var_out_start = tk.StringVar(value="1")
        self.var_out_end   = tk.StringVar(value="1")
        ttk.Entry(chan_frm, textvariable=self.var_out_start, width=4)\
            .grid(row=1, column=1)
        ttk.Label(chan_frm, text="–").grid(row=1, column=2)
        ttk.Entry(chan_frm, textvariable=self.var_out_end, width=4)\
            .grid(row=1, column=3)

        # Bottom buttons
        btn_frm = ttk.Frame(right)
        btn_frm.grid(row=4, column=0, sticky="ew", pady=6)
        self.btn_simulate = ttk.Button(btn_frm, text="Simulate",
                                       command=self._on_simulate)
        self.btn_simulate.pack(side=tk.LEFT, padx=4)
        self.btn_stop = ttk.Button(btn_frm, text="Stop",
                                   command=self._on_stop)
        self.btn_stop.pack(side=tk.LEFT, padx=4)
        ttk.Button(btn_frm, text="Edit",  command=self._on_edit)\
            .pack(side=tk.LEFT, padx=4)
        ttk.Button(btn_frm, text="Reset", command=self._on_reset)\
            .pack(side=tk.LEFT, padx=4)
        ttk.Button(btn_frm, text="Close", command=self._on_close)\
            .pack(side=tk.LEFT, padx=4)

    # ─────────────────────────────────────────────────────────────
    # Parameter display
    # ─────────────────────────────────────────────────────────────

    def refresh_params(self):
        """Update the parameter text display."""
        lines = [
            f"number of input channels  {g.nc}",
            f"center wavelength         {g.wl_c * 1e9:.2f} nm",
            f"channel spacing           {g.wl_d * 1e9:.2f} nm",
            f"waveguide width           {g.wg_w * 1e6:.2f} um",
            f"waveguide height          {g.wg_h * 1e6:.2f} um",
            f"polarization              {g.pol}",
            f"core index                {g.n_core:.6f}",
            f"clad index                {g.n_buff:.6f}",
            f"waveguide mode index      {g.n_eff0:.6f}",
            f"slab effective index      {g.n_eff1:.6f}",
            f"clad effective index      {g.n_eff2:.6f}",
            f"array waveguides          {g.na}",
            "",
            f"distance between adjacent array waveguides  {g.d_aw * 1e6:.6f} um",
            f"distance between adjacent i/o waveguides    {g.d_io * 1e6:.6f} um",
            "",
            f"phase order               {g.mphase:.0f}",
            f"path difference of awg    {g.d_path * 1e6:.6f} um",
            f"focal length              {g.f_l * 1e6:.4f} um",
            f"free spectral range       {g.fsr * 1e9:.4f} nm",
        ]
        text = "\n".join(lines)
        self.var_text.set(text)

        self.txt_params.configure(state="normal")
        self.txt_params.delete("1.0", tk.END)
        self.txt_params.insert("1.0", text)
        self.txt_params.configure(state="disabled")

    # ─────────────────────────────────────────────────────────────
    # Coupling coefficient buttons
    # ─────────────────────────────────────────────────────────────

    def _on_uniform(self):
        if g.na == 0:
            messagebox.showwarning("Warning", "Design not completed")
            return
        awg_calc.set_uniform_coupling()
        self.gauss_frm.grid_remove()
        self._show_coeff_plot()

    def _on_gauss(self):
        self.gauss_frm.grid()

    def _apply_gauss(self):
        if g.na == 0:
            messagebox.showwarning("Warning", "Design not completed")
            return
        try:
            width = float(self.var_gauss_w.get())
        except ValueError:
            width = 0.35
        awg_calc.set_gaussian_coupling(width)
        self._show_coeff_plot()

    def _on_bpm(self):
        messagebox.showinfo("BPM", "Load BPM coupling coefficients from a .res file.")
        filepath = filedialog.askopenfilename(
            title="BPM Result File",
            filetypes=[("BPM files", "*.res"), ("All Files", "*.*")])
        if not filepath:
            return
        self._load_bpm_coeffs(filepath)

    def _load_bpm_coeffs(self, filepath):
        """Load BPM amplitude & phase coupling coefficients from a .res file."""
        try:
            with open(filepath, "r") as fh:
                content = fh.read()

            # Find the marker line
            marker = "loss data in field amplitude and phase angle"
            idx = content.find(marker)
            if idx == -1:
                messagebox.showerror("Error", "BPM marker line not found in file")
                return

            # Skip 3 lines after marker, then read na amplitude + na phase values
            lines = content[idx + len(marker):].strip().splitlines()
            # skip 2 header lines + 1 data line
            data_start = 3
            values = []
            for line in lines[data_start:]:
                for tok in line.split():
                    values.append(float(tok))
                    if len(values) >= 2 * g.na:
                        break
                if len(values) >= 2 * g.na:
                    break

            for i in range(1, g.na + 1):
                g.Cpl_Coeff[i]  = values[i - 1] if i - 1 < len(values) else 0.0
            for i in range(1, g.na + 1):
                g.CplpCoeff[i] = values[g.na + i - 1] if g.na + i - 1 < len(values) else 0.0

            g.f_bpmcoeff = g.YES
            self._show_coeff_plot()

        except Exception as e:
            messagebox.showerror("BPM Load Error", str(e))

    def _show_coeff_plot(self):
        """Show coupling coefficient visualization."""
        from frm_cplcoeff import FrmCplCoeff
        if not hasattr(self, '_frm_cplcoeff') or self._frm_cplcoeff is None:
            self._frm_cplcoeff = FrmCplCoeff(self.master)
        self._frm_cplcoeff.show()
        self._frm_cplcoeff.plot_coeffs(g.Cpl_Coeff, g.na,
                                        gauss_width=float(self.var_gauss_w.get()))

    # ─────────────────────────────────────────────────────────────
    # Simulation
    # ─────────────────────────────────────────────────────────────

    def _on_simulate(self):
        if g.f_completed == g.NO_:
            messagebox.showwarning("Warning", "Design not completed")
            return
        if self._sim_running:
            return

        try:
            s_nw      = int(self.var_s_points.get())
            s_wls     = float(self.var_s_wlstart.get()) * 1e-9
            s_wle     = float(self.var_s_wlend.get())   * 1e-9
            s_wgd_w   = float(self.var_s_width.get())   * 1e-6
            s_wgd_h   = float(self.var_s_height.get())  * 1e-6
            i_in      = int(self.var_input_port.get())
            out_start = int(self.var_out_start.get())
            out_end   = int(self.var_out_end.get())
        except ValueError as e:
            messagebox.showerror("Input Error", str(e))
            return

        if s_nw < 10:
            messagebox.showwarning("Warning", "Too few points (min 10)")
            return

        output_ports = list(range(out_start, out_end + 1))

        # Set g.chk_1N_value so awg_calc knows the configuration
        g.chk_1N_value = (self.frm_variables.var_1n.get() == 1
                          if self.frm_variables else True)

        self._sim_stop   = False
        self._sim_running = True
        self.btn_simulate.configure(state="disabled")

        # Open simulation result window
        from frm_simulation import FrmSimulation
        if not hasattr(self, '_frm_sim') or self._frm_sim is None:
            self._frm_sim = FrmSimulation(self.master)
        self._frm_sim.prepare(s_wls * 1e9, s_wle * 1e9)
        self._frm_sim.show()

        def _run():
            try:
                results = awg_calc.simulate_spectral_response(
                    s_wls, s_wle, s_nw, s_wgd_w, s_wgd_h,
                    i_in, output_ports,
                    use_bpm=(g.f_bpmcoeff == g.YES),
                    progress_cb=self._sim_progress_cb,
                )
                # Save to xyz.dat
                with open("xyz.dat", "w") as fh:
                    for port, data in results.items():
                        for wl_nm, pwr in data:
                            fh.write(f"{wl_nm:.4f}  {pwr:.4f}\n")

                # Draw on simulation window
                self.master.after(0, lambda: self._frm_sim.draw_results(results))
            except Exception as e:
                self.master.after(0, lambda: messagebox.showerror("Simulation Error", str(e)))
            finally:
                self._sim_running = False
                self.master.after(0, lambda: self.btn_simulate.configure(state="normal"))

        t = threading.Thread(target=_run, daemon=True)
        t.start()

    def _sim_progress_cb(self, port, wl_nm, pwr_db):
        if self._sim_stop:
            raise InterruptedError("Simulation stopped")
        # Update wavelength display on sim window
        if hasattr(self, '_frm_sim') and self._frm_sim is not None:
            self.master.after(0, lambda wl=wl_nm:
                              self._frm_sim.update_wavelength(wl))

    def _on_stop(self):
        self._sim_stop = True

    def _on_reset(self):
        """Reset simulation variables to design defaults."""
        if self.frm_variables:
            self.var_s_width.set(self.frm_variables.var_w.get())
            self.var_s_height.set(self.frm_variables.var_h.get())
        self.var_s_loss.set("0")
        self.var_s_tempr.set("25")
        self.var_s_wlstart.set(f"{(g.wl_s - g.wl_d) * 1e9:.3f}")
        self.var_s_wlend.set(f"{(g.wl_e + g.wl_d) * 1e9:.3f}")
        self.var_s_points.set("500")
        # Update output port range
        self.var_out_start.set("1")
        self.var_out_end.set(str(g.nc))
        self.var_input_port.set("1")

    def _on_edit(self):
        if g.f_completed == g.NO_:
            messagebox.showwarning("Warning", "Design not completed")
            return
        from frm_coeff_edit import FrmCoeffEdit
        if not hasattr(self, '_frm_edit') or self._frm_edit is None:
            self._frm_edit = FrmCoeffEdit(self.master)
        self._frm_edit.load_data()
        self._frm_edit.show()

    def _on_close(self):
        self._sim_stop = True
        if self.window:
            self.window.withdraw()

    def _on_print(self):
        """Print parameters to console (no actual printer in Python port)."""
        print("=" * 60)
        print(g.VER_STRING)
        print(self.var_text.get())
        print("=" * 60)
