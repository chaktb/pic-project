"""
Refractive Index Calculation Form
Ported from frm_Index.frm (VB 5.0/6.0) to Python/Tkinter
"""

import os
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import math
import globals as g
import calc_slab as cs

# Directory that holds *.sel Sellmeier files
_DATA_DIR = os.path.join(os.path.dirname(__file__), "Data", "Sellmeier_C")


class FrmIndex:
    """Refractive Index Calculator window."""

    def __init__(self, master, frm_variables=None):
        self.master = master
        self.frm_variables = frm_variables
        self.window = None

        # Cached list of (display_name, filepath, data_dict) from Data/Sellmeier_C/
        self._sel_catalog = []
        # Status labels for each layer (type + WL range info)
        self._type_vars   = [tk.StringVar(value="") for _ in range(3)]

    def show(self):
        if self.window is None or not self.window.winfo_exists():
            self._build()
        self.window.deiconify()
        self.window.lift()

    def _build(self):
        # Scan .sel files every time the window is opened so new files are detected
        self._sel_catalog = cs.scan_sel_files(_DATA_DIR)

        self.window = tk.Toplevel(self.master)
        self.window.title("Refractive Index Calculation")
        self.window.resizable(False, False)
        self.window.protocol("WM_DELETE_WINDOW", self._on_cancel)

        self._build_menu()
        self._build_ui()
        self._load_tmp_data()

    def _build_menu(self):
        menubar = tk.Menu(self.window)
        file_m = tk.Menu(menubar, tearoff=0)
        file_m.add_command(label="New",     command=self._on_new)
        file_m.add_command(label="Open",    command=self._on_open)
        file_m.add_command(label="Save As", command=self._on_save_as)
        menubar.add_cascade(label="File", menu=file_m)

        wg_m = tk.Menu(menubar, tearoff=0)
        wg_m.add_command(label="Layer Data",              command=self._on_layer_data)
        wg_m.add_command(label="Read Mode Data",          command=self._on_read_mode_data)
        wg_m.add_command(label="Mode Calculation Batch",  command=self._on_mode_batch)
        menubar.add_cascade(label="Waveguide", menu=wg_m)

        self.window.config(menu=menubar)

    def _build_ui(self):
        # Combobox choice list from .sel catalog
        sel_names = ["(manual)"] + [name for name, _, _ in self._sel_catalog]

        # ── Top row: wavelength, temperature, polarization ──
        top = ttk.Frame(self.window, padding=6)
        top.grid(row=0, column=0, columnspan=3, sticky="ew")

        wl_frm = ttk.LabelFrame(top, text="wavelength", padding=4)
        wl_frm.grid(row=0, column=0, padx=4)
        self.var_wavelength = tk.StringVar()
        wl_entry = ttk.Entry(wl_frm, textvariable=self.var_wavelength, width=9)
        wl_entry.grid(row=0, column=0)
        ttk.Label(wl_frm, text="nm").grid(row=0, column=1)
        # Auto-recalc when wavelength changes
        self.var_wavelength.trace_add("write",
            lambda *_: self._auto_recalc_all())

        t_frm = ttk.LabelFrame(top, text="temperature", padding=4)
        t_frm.grid(row=0, column=1, padx=4)
        self.var_temp = tk.StringVar(value="25")
        ttk.Entry(t_frm, textvariable=self.var_temp, width=7).grid(row=0, column=0)
        ttk.Label(t_frm, text="°C").grid(row=0, column=1)
        self.var_temp.trace_add("write",
            lambda *_: self._auto_recalc_all())

        pol_frm = ttk.LabelFrame(top, text="Polarization", padding=4)
        pol_frm.grid(row=0, column=2, padx=4)
        self.var_pol = tk.StringVar(value="X")
        ttk.Combobox(pol_frm, textvariable=self.var_pol,
                     values=["X", "Y", "S"], width=4).grid(row=0, column=0)

        # ── Material selection + index display (one frame per layer) ──
        mat_frame = ttk.Frame(self.window, padding=6)
        mat_frame.grid(row=1, column=0, columnspan=3, sticky="ew")

        self.var_mtl       = [tk.StringVar() for _ in range(3)]
        self.var_indx      = [tk.StringVar() for _ in range(3)]
        self.var_sel_mat   = [tk.StringVar(value="(manual)") for _ in range(3)]
        # store path of currently loaded .sel per layer (None = manual)
        self._sel_loaded   = [None, None, None]

        layer_labels = [("core material", 0),
                        ("over material", 1),
                        ("buffer material", 2)]

        for col, (label, layer) in enumerate(layer_labels):
            f = ttk.LabelFrame(mat_frame, text=label, padding=6)
            f.grid(row=0, column=col, padx=4, pady=4, sticky="nsew")

            # Row 0: material selector combobox
            ttk.Label(f, text="Material:").grid(row=0, column=0, sticky="w")
            cmb = ttk.Combobox(f, textvariable=self.var_sel_mat[layer],
                               values=sel_names, state="readonly", width=16)
            cmb.grid(row=0, column=1, columnspan=2, sticky="ew", pady=2)
            cmb.bind("<<ComboboxSelected>>",
                     lambda e, lay=layer: self._on_sel_mat_change(lay))

            # Row 1: manual name entry (shown when "(manual)" selected)
            ttk.Label(f, text="Name:").grid(row=1, column=0, sticky="w")
            ttk.Entry(f, textvariable=self.var_mtl[layer], width=16)\
                .grid(row=1, column=1, columnspan=2, sticky="ew")

            # Row 2: calculated index (read-only display)
            ttk.Label(f, text="n:").grid(row=2, column=0, sticky="w")
            ttk.Entry(f, textvariable=self.var_indx[layer],
                      width=16, state="readonly")\
                .grid(row=2, column=1, columnspan=2, sticky="ew")

            # Row 3: type label
            self._build_type_label(f, layer, row=3)

        # ── Waveguide dimensions ──
        dim_frame = ttk.Frame(self.window, padding=6)
        dim_frame.grid(row=2, column=0, columnspan=2, sticky="ew")

        h_frm = ttk.LabelFrame(dim_frame, text="core height", padding=4)
        h_frm.grid(row=0, column=0, padx=4)
        self.var_hcore = tk.StringVar()
        ttk.Entry(h_frm, textvariable=self.var_hcore, width=8).grid(row=0, column=0)
        ttk.Label(h_frm, text="μm").grid(row=0, column=1)

        w_frm = ttk.LabelFrame(dim_frame, text="core width", padding=4)
        w_frm.grid(row=0, column=1, padx=4)
        self.var_wcore = tk.StringVar()
        ttk.Entry(w_frm, textvariable=self.var_wcore, width=8).grid(row=0, column=0)
        ttk.Label(w_frm, text="μm").grid(row=0, column=1)

        hbuf_frm = ttk.LabelFrame(dim_frame, text="buffer height", padding=4)
        hbuf_frm.grid(row=0, column=2, padx=4)
        self.var_hbuff = tk.StringVar()
        ttk.Entry(hbuf_frm, textvariable=self.var_hbuff, width=8).grid(row=0, column=0)
        ttk.Label(hbuf_frm, text="μm").grid(row=0, column=1)

        # ── Effective index results + EIM method ──
        eff_frm = ttk.LabelFrame(self.window, text="eff. indx. Re", padding=4)
        eff_frm.grid(row=2, column=2, padx=6, pady=4, sticky="nsew")

        self.var_ercore = tk.StringVar()
        self.var_erclad = tk.StringVar()
        self.var_betta  = tk.StringVar()

        ttk.Label(eff_frm, text="Core").grid(row=0, column=0, sticky="e")
        ttk.Entry(eff_frm, textvariable=self.var_ercore, width=14, state="readonly")\
            .grid(row=0, column=1)
        ttk.Label(eff_frm, text="Clad").grid(row=1, column=0, sticky="e")
        ttk.Entry(eff_frm, textvariable=self.var_erclad, width=14, state="readonly")\
            .grid(row=1, column=1)
        ttk.Label(eff_frm, text="Beta").grid(row=2, column=0, sticky="e")
        ttk.Entry(eff_frm, textvariable=self.var_betta,  width=14, state="readonly")\
            .grid(row=2, column=1)

        self.var_eim = tk.IntVar(value=0)
        for i, txt in enumerate(["EIM", "IEIM", "HEIM"]):
            ttk.Radiobutton(eff_frm, text=txt, variable=self.var_eim, value=i)\
                .grid(row=3, column=i, padx=2)

        # ── Mode index table ──
        table_frm = ttk.LabelFrame(self.window, text="Mode Index Table", padding=4)
        table_frm.grid(row=0, column=3, rowspan=3, padx=6, pady=4, sticky="nsew")

        self.mode_table = ttk.Treeview(table_frm, columns=("no","wl","idx"),
                                        show="headings", height=8)
        self.mode_table.heading("no",  text="NO")
        self.mode_table.heading("wl",  text="wave L (nm)")
        self.mode_table.heading("idx", text="index")
        self.mode_table.column("no",  width=35)
        self.mode_table.column("wl",  width=80)
        self.mode_table.column("idx", width=100)
        self.mode_table.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        scroll = ttk.Scrollbar(table_frm, orient=tk.VERTICAL,
                               command=self.mode_table.yview)
        scroll.pack(side=tk.RIGHT, fill=tk.Y)
        self.mode_table.configure(yscrollcommand=scroll.set)

        rng_frm = ttk.Frame(table_frm)
        rng_frm.pack(fill=tk.X, pady=4)
        ttk.Label(rng_frm, text="WL Range (nm):").grid(row=0, column=0, sticky="w")
        self.var_wlmin  = tk.StringVar()
        self.var_wlmax  = tk.StringVar()
        self.var_nwavel = tk.StringVar(value="10")
        ttk.Entry(rng_frm, textvariable=self.var_wlmin,  width=6).grid(row=0, column=1)
        ttk.Label(rng_frm, text="–").grid(row=0, column=2)
        ttk.Entry(rng_frm, textvariable=self.var_wlmax,  width=6).grid(row=0, column=3)
        ttk.Entry(rng_frm, textvariable=self.var_nwavel, width=4).grid(row=0, column=4)

        self.var_tcoeff = tk.StringVar(value="0.0")
        ttk.Label(rng_frm, text="Temp. Coeff:").grid(row=1, column=0, sticky="w")
        ttk.Entry(rng_frm, textvariable=self.var_tcoeff, width=10)\
            .grid(row=1, column=1, columnspan=3, sticky="w")
        ttk.Label(rng_frm, text="/°C").grid(row=1, column=4)

        # ── Action buttons ──
        btn_frm = ttk.Frame(self.window, padding=6)
        btn_frm.grid(row=3, column=0, columnspan=4, sticky="e")

        ttk.Button(btn_frm, text="Calculate",    command=self._on_calculate)\
            .pack(side=tk.LEFT, padx=4)
        ttk.Button(btn_frm, text="Mode Profile", command=self._on_mode_profile)\
            .pack(side=tk.LEFT, padx=4)
        ttk.Button(btn_frm, text="OK",           command=self._on_ok)\
            .pack(side=tk.LEFT, padx=4)
        ttk.Button(btn_frm, text="Cancel",       command=self._on_cancel)\
            .pack(side=tk.LEFT, padx=4)

    def _build_type_label(self, parent, layer, row):
        """Small status label showing type + WL range of loaded .sel."""
        lbl = ttk.Label(parent, textvariable=self._type_vars[layer],
                        foreground="#888888", font=("", 8))
        lbl.grid(row=row, column=0, columnspan=3, sticky="w")

    # ─────────────────────────────────────────────────────────────
    # .sel material selector
    # ─────────────────────────────────────────────────────────────

    def _on_sel_mat_change(self, layer):
        """Called when user picks a material from the Combobox."""
        choice = self.var_sel_mat[layer].get()
        if choice == "(manual)":
            self._sel_loaded[layer] = None
            self._type_vars[layer].set("")
            return

        # Find matching entry in catalog
        for name, filepath, data in self._sel_catalog:
            if name == choice:
                self._apply_sel(layer, data, filepath)
                return

    def _apply_sel(self, layer, data, filepath):
        """Load .sel data into globals and auto-calculate index."""
        cs.apply_sel_to_globals(layer, data)
        self._sel_loaded[layer] = filepath

        # Update material name display
        self.var_mtl[layer].set(data["name"])

        # Update type/range status label
        typ = data["type"]
        wmin, wmax = data["wl_min"], data["wl_max"]
        if wmin or wmax:
            self._type_vars[layer].set(
                f"{typ}  {int(wmin)}–{int(wmax)} nm  Tc={data['tcoeff']:.2e}/°C")
        else:
            self._type_vars[layer].set(f"{typ}  Tc={data['tcoeff']:.2e}/°C")

        # Auto-calculate index at current wavelength
        self._recalc_layer_index(layer)

    def _recalc_layer_index(self, layer):
        """Recalculate n for one layer using Sellmeier and update var_indx."""
        try:
            wl_nm = float(self.var_wavelength.get())
            temp  = float(self.var_temp.get() or "25")
        except ValueError:
            return
        if wl_nm <= 0:
            return
        try:
            n = cs.calc_sellmeier(layer, wl_nm, temp)
            self.var_indx[layer].set(f"{n:.10f}")
        except Exception:
            pass

    def _auto_recalc_all(self):
        """Recalculate all layers that have a .sel loaded."""
        for layer in range(3):
            if self._sel_loaded[layer] is not None:
                self._recalc_layer_index(layer)

    # ─────────────────────────────────────────────────────────────
    # Actions
    # ─────────────────────────────────────────────────────────────

    def _on_calculate(self):
        try:
            wl_nm  = float(self.var_wavelength.get())
            pol    = self.var_pol.get().upper()
            n_core = float(self.var_indx[0].get())
            n_over = float(self.var_indx[1].get())
            n_buff = float(self.var_indx[2].get())
            wgd_h  = float(self.var_hcore.get()) * 1e-6
            wgd_w  = float(self.var_wcore.get()) * 1e-6
            buf_h  = float(self.var_hbuff.get()) * 1e-6 if self.var_hbuff.get() else 0.0
            wl     = wl_nm * 1e-9
        except ValueError as e:
            messagebox.showerror("Input Error", str(e))
            return

        if wl <= 0:
            messagebox.showerror("Error", "Wavelength not defined")
            return
        if n_core <= max(n_over, n_buff):
            messagebox.showerror("Error", "Core index too low")
            return

        method = self.var_eim.get()
        try:
            if method == 0:  # EIM
                n_mode, n_side, betta = cs.Calc_NEIM(
                    wl, pol, n_core, n_over, n_buff, wgd_h, wgd_w, buf_h,
                    None, None, None)
            elif method == 2:  # HEIM
                # Get n_mode from table
                n_init = cs.Calc_IndxTable(4, wl_nm)
                if n_init <= 0:
                    n_init = n_core
                temp = float(self.var_temp.get()) if self.var_temp.get() else 25.0
                tc   = float(self.var_tcoeff.get()) if self.var_tcoeff.get() else 0.0
                n_init += tc * (temp - 25.0)
                n_mode, n_side, betta = cs.Calc_HEIM(
                    wl, pol, n_core, n_over, n_buff, wgd_h, wgd_w, buf_h,
                    None, None, None, n_init)
            else:  # IEIM – not implemented
                messagebox.showinfo("Info", "IEIM method not implemented")
                return
        except Exception as e:
            messagebox.showerror("Calculation Error", str(e))
            return

        self.var_ercore.set(f"{n_mode:.10f}")
        self.var_erclad.set(f"{n_side:.10f}")
        self.var_betta.set(f"{betta:.4f}")

    def _on_mode_profile(self):
        try:
            wl_nm  = float(self.var_wavelength.get())
            pol    = self.var_pol.get().upper()
            wgd_w  = float(self.var_wcore.get()) * 1e-6
            n_core = float(self.var_ercore.get())
            n_side = float(self.var_erclad.get())
            wl     = wl_nm * 1e-9
        except ValueError as e:
            messagebox.showerror("Input Error", str(e))
            return

        filepath = filedialog.asksaveasfilename(
            title="Save Mode Profile",
            defaultextension=".mod",
            filetypes=[("Mode Profile", "*.mod"), ("All Files", "*.*")])
        if not filepath:
            return

        try:
            cs.Save_ModePrf(filepath, wl, pol, wgd_w, n_core, n_side)
            messagebox.showinfo("Done", f"Mode profile saved to {filepath}")
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def _on_ok(self):
        """Apply calculated values back to main form."""
        if self.frm_variables:
            try:
                n_core = float(self.var_indx[0].get())
                n_over = float(self.var_indx[1].get())
                n_buff = float(self.var_indx[2].get())
            except ValueError:
                n_core = n_over = n_buff = 0.0

            self.frm_variables.var_wlc.set(self.var_wavelength.get())
            self.frm_variables.var_w.set(self.var_wcore.get())
            self.frm_variables.var_h.set(self.var_hcore.get())
            self.frm_variables.var_pol.set(self.var_pol.get())
            self.frm_variables.update_index_display(n_core, n_over, n_buff)
        self.window.withdraw()

    def _on_cancel(self):
        self.window.withdraw()

    def _on_new(self):
        for v in self.var_mtl + self.var_indx:
            v.set("")
        self.var_ercore.set("")
        self.var_erclad.set("")
        self.var_betta.set("")

    def _on_open(self):
        filepath = filedialog.askopenfilename(
            title="Open Index File",
            filetypes=[("Index File", "*.inx"), ("All Files", "*.*")])
        if not filepath:
            return
        self._load_inx_file(filepath)

    def _on_save_as(self):
        filepath = filedialog.asksaveasfilename(
            title="Save Index File",
            defaultextension=".inx",
            filetypes=[("Index File", "*.inx"), ("All Files", "*.*")])
        if not filepath:
            return
        self._save_inx_file(filepath)

    def _on_layer_data(self):
        from frm_rindex import FrmRIndex
        if not hasattr(self, "_frm_rindex") or self._frm_rindex is None:
            self._frm_rindex = FrmRIndex(self.master, frm_index=self)
        self._frm_rindex.show()

    def _on_read_mode_data(self):
        filepath = filedialog.askopenfilename(
            title="Open Mode Index Table",
            filetypes=[("Mode Index Table", "*.int"), ("All Files", "*.*")])
        if not filepath:
            return
        self._load_mode_table(filepath)

    def _on_mode_batch(self):
        messagebox.showinfo("Info", "Mode Calculation Batch: use menu Waveguide > Layer Data to set materials, then use Mode Index Table to generate data.")

    def _load_mode_table(self, filepath):
        """Load mode index table from .int file."""
        try:
            rows = []
            with open(filepath, "r") as fh:
                for line in fh:
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        rows.append((float(parts[0]), float(parts[1])))
            # Update global table (layer 4 = mode table)
            g.nr_table[4][0] = [r[0] for r in rows]
            g.nr_table[4][1] = [r[1] for r in rows]
            g.nr_index[4] = len(rows)
            g.nr_Tcoeff[4] = float(self.var_tcoeff.get()) if self.var_tcoeff.get() else 0.0

            # Update treeview
            for item in self.mode_table.get_children():
                self.mode_table.delete(item)
            for i, (wl, idx) in enumerate(rows):
                self.mode_table.insert("", "end", values=(i+1, f"{wl:.3f}", f"{idx:.8f}"))
        except Exception as e:
            messagebox.showerror("Load Error", str(e))

    def _load_inx_file(self, filepath):
        """Load .inx refractive index file."""
        try:
            with open(filepath, "r") as fh:
                lines = fh.readlines()
            if len(lines) >= 3:
                for i in range(3):
                    parts = lines[i].strip().split(",")
                    self.var_mtl[i].set(parts[0] if len(parts) > 0 else "")
                    self.var_indx[i].set(parts[1] if len(parts) > 1 else "")
        except Exception as e:
            messagebox.showerror("Load Error", str(e))

    def _save_inx_file(self, filepath):
        try:
            with open(filepath, "w") as fh:
                for i in range(3):
                    fh.write(f"{self.var_mtl[i].get()},{self.var_indx[i].get()}\n")
        except Exception as e:
            messagebox.showerror("Save Error", str(e))

    def _load_tmp_data(self):
        """Restore previously entered data."""
        try:
            self.var_wavelength.set(str(g.wl_c * 1e9) if g.wl_c else "")
            self.var_hcore.set(str(g.wg_h * 1e6) if g.wg_h else "")
            self.var_wcore.set(str(g.wg_w * 1e6) if g.wg_w else "")
            self.var_pol.set(g.pol or "X")
            self.var_indx[0].set(str(g.n_core) if g.n_core else "")
            self.var_indx[1].set(str(g.n_over) if g.n_over else "")
            self.var_indx[2].set(str(g.n_buff) if g.n_buff else "")
        except Exception:
            pass

    def set_inputs(self, wl_nm="", w_um="", h_um="", pol="X",
                   n_core="", n_over="", n_buff=""):
        """Set input fields from main form."""
        if self.window is None or not self.window.winfo_exists():
            self._build()
        self.var_wavelength.set(wl_nm)
        self.var_wcore.set(w_um)
        self.var_hcore.set(h_um)
        self.var_pol.set(pol)
        self.var_indx[0].set(n_core)
        self.var_indx[1].set(n_over)
        self.var_indx[2].set(n_buff)

    def mode_calc(self, wl_m):
        """
        Calculate mode indices at wl_m [m] using current window settings.
        Returns (n_mode, n_slab, n_side).
        Called from awg_calc.calc_fsr().
        """
        import awg_calc
        try:
            n_core = float(self.var_indx[0].get())
            n_over = float(self.var_indx[1].get())
            n_buff = float(self.var_indx[2].get())
            wgd_h  = float(self.var_hcore.get()) * 1e-6
            wgd_w  = float(self.var_wcore.get()) * 1e-6
            pol    = self.var_pol.get().upper() or "X"
        except ValueError:
            return awg_calc.mode_calc(wl_m)

        return awg_calc.mode_calc(
            wl_m,
            n_core_val=n_core,
            n_over_val=n_over,
            n_buff_val=n_buff,
            wgd_w=wgd_w,
            wgd_h=wgd_h,
            pol_val=pol,
        )
