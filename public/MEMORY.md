"""
Refractive Index Data Management Form
Ported from frm_RIndex.frm (VB 5.0/6.0) to Python/Tkinter

Manages Sellmeier equation coefficients and temperature dependence
for core, overcladding, and buffer materials.
Supports .sel material files and index delta (%) correction.
"""

import os
import tkinter as tk
from tkinter import ttk, messagebox
import math
import globals as g
import calc_slab as cs

_DATA_DIR = os.path.join(os.path.dirname(__file__), "Data", "Sellmeier_C")


class FrmRIndex:
    """Refractive index model editor."""

    def __init__(self, master, frm_index=None):
        self.master = master
        self.frm_index = frm_index   # reference to FrmIndex for propagation
        self.window = None

        # .sel catalog: list of (display_name, filepath, data_dict)
        self._sel_catalog = []

        # Material names
        self.var_mtl_name  = [tk.StringVar() for _ in range(3)]
        # Refractive index output (calculated)
        self.var_rindex    = [tk.StringVar() for _ in range(3)]
        # Index delta (%)
        self.var_delta     = [tk.StringVar(value="0.0") for _ in range(3)]
        # Wavelength & temperature
        self.var_wavelength = tk.StringVar(value="1550")
        self.var_temp       = tk.StringVar(value="25")

        # Sellmeier type per material
        self.var_type = [tk.StringVar(value="B") for _ in range(3)]
        # Sellmeier coefficients [material][coeff_index]
        self.var_coeffs = [[tk.StringVar(value="0.0") for _ in range(6)]
                           for _ in range(3)]
        # Temperature coefficient
        self.var_tcoeff = [tk.StringVar(value="0.0") for _ in range(3)]

        # Selected material from combobox per layer
        self.var_sel_mat = [tk.StringVar(value="(manual)") for _ in range(3)]
        # Loaded filepath per layer (None = manual)
        self._sel_loaded = [None, None, None]

    def show(self):
        if self.window is None or not self.window.winfo_exists():
            self._sel_catalog = cs.scan_sel_files(_DATA_DIR)
            self._build()
        self.window.deiconify()
        self.window.lift()
        # Sync wavelength from frm_index if available
        self._sync_wavelength_from_index()

    def _sync_wavelength_from_index(self):
        """Pull current wavelength/temp from frm_index if it exists."""
        if self.frm_index is None:
            return
        fi = self.frm_index
        if hasattr(fi, "var_wavelength"):
            try:
                wl = fi.var_wavelength.get()
                if wl:
                    self.var_wavelength.set(wl)
            except Exception:
                pass
        if hasattr(fi, "var_temp"):
            try:
                t = fi.var_temp.get()
                if t:
                    self.var_temp.set(t)
            except Exception:
                pass

    def _build(self):
        self.window = tk.Toplevel(self.master)
        self.window.title("Refractive Index Data")
        self.window.resizable(False, False)
        self.window.protocol("WM_DELETE_WINDOW", self._on_close)
        self._build_ui()
        # Set up auto-recalculation traces
        self.var_wavelength.trace_add("write", lambda *_: self._auto_recalc_all())
        self.var_temp.trace_add("write",       lambda *_: self._auto_recalc_all())

    def _build_ui(self):
        sel_names = ["(manual)"] + [name for name, _, _ in self._sel_catalog]

        # ── Top row: Wavelength & Temperature ──
        top = ttk.Frame(self.window, padding=6)
        top.grid(row=0, column=0, columnspan=3, sticky="ew")

        ttk.Label(top, text="Wavelength (nm):").grid(row=0, column=0, sticky="w")
        ttk.Entry(top, textvariable=self.var_wavelength, width=10)\
            .grid(row=0, column=1)
        ttk.Label(top, text="Temperature (°C):").grid(row=0, column=2, sticky="w",
                                                       padx=(16, 0))
        ttk.Entry(top, textvariable=self.var_temp, width=8)\
            .grid(row=0, column=3)

        # ── Material frames ──
        mat_labels = ["Core", "Overcladding", "Buffer"]
        for col, label in enumerate(mat_labels):
            self._build_material_frame(col, label, sel_names)

        # ── Refractive index results ──
        res_frm = ttk.Frame(self.window, padding=6)
        res_frm.grid(row=5, column=0, columnspan=3, sticky="ew")

        for col, (name, var) in enumerate(zip(mat_labels, self.var_rindex)):
            ttk.Label(res_frm, text=f"n({name}):").grid(
                row=0, column=col * 2, sticky="w", padx=4)
            ttk.Entry(res_frm, textvariable=var, width=14, state="readonly")\
                .grid(row=0, column=col * 2 + 1, sticky="w")

        # ── Buttons ──
        btn_frm = ttk.Frame(self.window, padding=6)
        btn_frm.grid(row=6, column=0, columnspan=3, sticky="e")

        ttk.Button(btn_frm, text="Calculate All",
                   command=self.cmd_CalcAll_Click)\
            .pack(side=tk.LEFT, padx=4)
        ttk.Button(btn_frm, text="Apply to Index Calc",
                   command=self._on_apply_to_index)\
            .pack(side=tk.LEFT, padx=4)
        ttk.Button(btn_frm, text="OK",
                   command=self.cmd_OK_Click)\
            .pack(side=tk.LEFT, padx=4)
        ttk.Button(btn_frm, text="Cancel",
                   command=self._on_close)\
            .pack(side=tk.LEFT, padx=4)

    def _build_material_frame(self, col, label, sel_names):
        frm = ttk.LabelFrame(self.window, text=label, padding=6)
        frm.grid(row=1, column=col, padx=6, pady=4, sticky="nsew")

        # Row 0: material selector combobox
        ttk.Label(frm, text="Material:").grid(row=0, column=0, sticky="w")
        cmb = ttk.Combobox(frm, textvariable=self.var_sel_mat[col],
                           values=sel_names, state="readonly", width=16)
        cmb.grid(row=0, column=1, columnspan=2, sticky="ew", pady=2)
        cmb.bind("<<ComboboxSelected>>",
                 lambda e, c=col: self._on_sel_mat_change(c))

        # Row 1: material name (auto-filled or manual)
        ttk.Label(frm, text="Name:").grid(row=1, column=0, sticky="w")
        ttk.Entry(frm, textvariable=self.var_mtl_name[col], width=16)\
            .grid(row=1, column=1, columnspan=2, sticky="ew")

        # Row 2: type combobox
        ttk.Label(frm, text="Type:").grid(row=2, column=0, sticky="w")
        ttk.Combobox(frm, textvariable=self.var_type[col],
                     values=["A", "B"], width=5)\
            .grid(row=2, column=1, sticky="w")

        # Rows 3–8: Sellmeier coefficients c1..c6
        for ci in range(6):
            ttk.Label(frm, text=f"c{ci+1}:").grid(row=3 + ci, column=0, sticky="w")
            ttk.Entry(frm, textvariable=self.var_coeffs[col][ci], width=14)\
                .grid(row=3 + ci, column=1, columnspan=2)

        # Row 9: temperature coefficient
        ttk.Label(frm, text="Temp.coeff (/°C):").grid(
            row=9, column=0, sticky="w", pady=(4, 0))
        ttk.Entry(frm, textvariable=self.var_tcoeff[col], width=12)\
            .grid(row=9, column=1, columnspan=2)

        # Row 10: index delta (%)
        ttk.Label(frm, text="Δn (%):").grid(row=10, column=0, sticky="w", pady=(4, 0))
        delta_e = ttk.Entry(frm, textvariable=self.var_delta[col], width=10)
        delta_e.grid(row=10, column=1, columnspan=2, sticky="w")
        # Re-calculate when delta changes
        self.var_delta[col].trace_add("write",
            lambda *_, c=col: self._recalc_index(c))

    # ─────────────────────────────────────────────────────────────
    # .sel material selector
    # ─────────────────────────────────────────────────────────────

    def _on_sel_mat_change(self, col):
        choice = self.var_sel_mat[col].get()
        if choice == "(manual)":
            self._sel_loaded[col] = None
            return

        for name, filepath, data in self._sel_catalog:
            if name == choice:
                self._apply_sel(col, data, filepath)
                return

    def _apply_sel(self, col, data, filepath):
        """Load .sel coefficients into this form and recalculate."""
        cs.apply_sel_to_globals(col, data)
        self._sel_loaded[col] = filepath

        self.var_mtl_name[col].set(data["name"])
        self.var_type[col].set(data["type"])
        self.var_tcoeff[col].set(str(data["tcoeff"]))

        # Fill coefficient fields
        coeffs = data["coeffs"]
        for ci in range(6):
            val = coeffs[ci] if ci < len(coeffs) else 0.0
            self.var_coeffs[col][ci].set(str(val))

        self._recalc_index(col)

    # ─────────────────────────────────────────────────────────────
    # Index calculation
    # ─────────────────────────────────────────────────────────────

    def _recalc_index(self, col):
        """Recalculate n for one material column and update var_rindex."""
        try:
            wl_nm = float(self.var_wavelength.get())
            temp  = float(self.var_temp.get() or "25")
        except ValueError:
            return
        if wl_nm <= 0:
            return

        n = self._calc_index(col, wl_nm, temp)
        self.var_rindex[col].set(f"{n:.8f}")

    def _auto_recalc_all(self):
        """Recalculate all columns that have a .sel loaded or have coefficients."""
        for col in range(3):
            # Only auto-recalc if a .sel is loaded; manual entry requires explicit button
            if self._sel_loaded[col] is not None:
                self._recalc_index(col)

    def cmd_CalcAll_Click(self):
        """Calculate refractive indices for all materials."""
        try:
            wl_nm = float(self.var_wavelength.get())
            temp  = float(self.var_temp.get())
        except ValueError:
            messagebox.showerror("Error", "Invalid wavelength or temperature")
            return

        for col in range(3):
            n = self._calc_index(col, wl_nm, temp)
            self.var_rindex[col].set(f"{n:.8f}")
            # Store Sellmeier data in globals
            g.sellmeier[col]["name"]   = self.var_mtl_name[col].get()
            g.sellmeier[col]["type"]   = self.var_type[col].get()
            g.sellmeier[col]["coeffs"] = [
                self._flt(self.var_coeffs[col][i]) for i in range(6)]
            g.nr_Tcoeff[col] = self._flt(self.var_tcoeff[col])

    def _calc_index(self, col, wl_nm, temp):
        """Calculate refractive index for material col at wl_nm, temp,
        then apply delta correction: n_final = n / (1 - delta/100)."""
        type_ = self.var_type[col].get()
        c = [self._flt(self.var_coeffs[col][i]) for i in range(6)]
        tc = self._flt(self.var_tcoeff[col])
        wl_um = wl_nm / 1000.0

        if type_ in ("A", "TypeA"):
            # TypeA: n² = 1 + Σ Bi·λ²/(λ²−Li²)  (Li squared)
            n2 = 1.0
            for k in range(3):
                B = c[k * 2]
                L = c[k * 2 + 1]
                if L != 0:
                    n2 += B * wl_um ** 2 / (wl_um ** 2 - L ** 2)
            n = math.sqrt(max(n2, 1.0))
        elif type_ in ("B", "TypeB"):
            # TypeB: n² = 1 + Σ Bi·λ²/(λ²−Li)   (Li NOT squared)
            n2 = 1.0
            for k in range(3):
                B = c[k * 2]
                L = c[k * 2 + 1]
                if L != 0:
                    n2 += B * wl_um ** 2 / (wl_um ** 2 - L)
            n = math.sqrt(max(n2, 1.0))
        else:
            n = c[0] if c[0] > 0 else 1.44

        # Temperature correction
        n += tc * (temp - 25.0)

        # Index delta correction: n_final = n / (1 - delta/100)
        try:
            delta = float(self.var_delta[col].get())
        except (ValueError, AttributeError):
            delta = 0.0
        if delta != 0.0:
            denom = 1.0 - delta / 100.0
            if abs(denom) > 1e-10:
                n = n / denom

        return n

    def _flt(self, var):
        try:
            return float(var.get())
        except (ValueError, AttributeError):
            return 0.0

    # ─────────────────────────────────────────────────────────────
    # Actions
    # ─────────────────────────────────────────────────────────────

    def _on_apply_to_index(self):
        """Push calculated indices to the open frm_index window."""
        if self.frm_index is None:
            messagebox.showinfo("Info",
                "No Refractive Index Calculation window linked.")
            return

        # Ensure indices are calculated
        try:
            wl_nm = float(self.var_wavelength.get())
            temp  = float(self.var_temp.get())
        except ValueError:
            messagebox.showerror("Error", "Invalid wavelength or temperature")
            return

        for col in range(3):
            n = self._calc_index(col, wl_nm, temp)
            self.var_rindex[col].set(f"{n:.8f}")

        fi = self.frm_index
        # Make sure window is built
        if fi.window is None or not fi.window.winfo_exists():
            fi.show()

        # Propagate material names and indices
        for col in range(3):
            fi.var_mtl[col].set(self.var_mtl_name[col].get())
            fi.var_indx[col].set(self.var_rindex[col].get())

        # Propagate Sellmeier to globals
        g.n_core = self._calc_index(0, wl_nm, temp)
        g.n_over = self._calc_index(1, wl_nm, temp)
        g.n_buff = self._calc_index(2, wl_nm, temp)

        messagebox.showinfo("Applied",
            f"n_core={g.n_core:.8f}\n"
            f"n_over={g.n_over:.8f}\n"
            f"n_buff={g.n_buff:.8f}\n\nApplied to Refractive Index Calculation.")

    def cmd_OK_Click(self):
        """Apply calculated indices to globals and propagate."""
        try:
            wl_nm = float(self.var_wavelength.get())
            temp  = float(self.var_temp.get())
        except ValueError:
            wl_nm, temp = 1550.0, 25.0

        n_core = self._calc_index(0, wl_nm, temp)
        n_over = self._calc_index(1, wl_nm, temp)
        n_buff = self._calc_index(2, wl_nm, temp)

        g.n_core = n_core
        g.n_over = n_over
        g.n_buff = n_buff

        # Update Sellmeier globals
        for col in range(3):
            g.sellmeier[col]["name"]   = self.var_mtl_name[col].get()
            g.sellmeier[col]["type"]   = self.var_type[col].get()
            g.sellmeier[col]["coeffs"] = [
                self._flt(self.var_coeffs[col][i]) for i in range(6)]
            g.nr_Tcoeff[col] = self._flt(self.var_tcoeff[col])

        # Propagate to frm_index if open
        if self.frm_index is not None:
            fi = self.frm_index
            if fi.window and fi.window.winfo_exists():
                for col in range(3):
                    fi.var_mtl[col].set(self.var_mtl_name[col].get())
                    fi.var_indx[col].set(self.var_rindex[col].get())

        self._on_close()

    def Save_Tmp_Data(self):
        """Store temp data from frm_index."""
        pass

    def Recover_Tmp_Data(self):
        """Restore temp data."""
        pass

    def _on_close(self):
        if self.window:
            self.window.withdraw()

    @property
    def txt_RIndex(self):
        """Compatibility shim: txt_RIndex[i].Text → var_rindex[i]."""
        class _Entry:
            def __init__(self, var):
                self._var = var

            @property
            def Text(self):
                return self._var.get()

        return [_Entry(self.var_rindex[i]) for i in range(3)]
