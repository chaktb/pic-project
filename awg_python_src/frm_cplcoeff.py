"""
AWG Loss & Parameter Editor Form
Ported from frm_CoeffEdit.frm (VB 5.0/6.0) to Python/Tkinter

Grid-based editing of per-waveguide AWG parameters:
Loss, Phase, Length, Bend radius, Taper length, Width, Height
"""

import tkinter as tk
from tkinter import ttk, messagebox
import math
import globals as g


class FrmCoeffEdit:
    """AWG waveguide parameter editor."""

    COL_HEADERS = ["#", "Loss (dB/cm)", "Phase (deg)", "Length (um)",
                   "Arc R (um)", "Arc L (um)", "Width (um)", "Height (um)"]

    def __init__(self, master):
        self.master = master
        self.window = None

    def show(self):
        if self.window is None or not self.window.winfo_exists():
            self._build()
        self.window.deiconify()
        self.window.lift()

    def _build(self):
        self.window = tk.Toplevel(self.master)
        self.window.title("AWG Waveguide Parameters Editor")
        self.window.resizable(True, True)
        self.window.protocol("WM_DELETE_WINDOW", self._on_close)
        self._build_ui()

    def _build_ui(self):
        # Bulk-set controls
        ctrl = ttk.Frame(self.window, padding=6)
        ctrl.pack(fill=tk.X)

        ttk.Label(ctrl, text="Loss (dB/cm):").pack(side=tk.LEFT)
        self.var_loss = tk.StringVar(value="0")
        ttk.Entry(ctrl, textvariable=self.var_loss, width=8)\
            .pack(side=tk.LEFT, padx=2)

        ttk.Button(ctrl, text="Apply Loss",   command=self._apply_loss)\
            .pack(side=tk.LEFT, padx=4)
        ttk.Button(ctrl, text="Zero Phase",   command=self._zero_phase)\
            .pack(side=tk.LEFT, padx=4)
        ttk.Button(ctrl, text="Random Phase", command=self._random_phase)\
            .pack(side=tk.LEFT, padx=4)

        # Grid (Treeview)
        tree_frm = ttk.Frame(self.window)
        tree_frm.pack(fill=tk.BOTH, expand=True, padx=6, pady=6)

        self.tree = ttk.Treeview(
            tree_frm,
            columns=["no", "loss", "phase", "length",
                     "arc_r", "arc_l", "width", "height"],
            show="headings",
            height=20,
        )
        for col, header in zip(
                ["no", "loss", "phase", "length", "arc_r", "arc_l", "width", "height"],
                self.COL_HEADERS):
            self.tree.heading(col, text=header)
            self.tree.column(col, width=80, anchor="center")
        self.tree.column("no", width=40)

        vsb = ttk.Scrollbar(tree_frm, orient=tk.VERTICAL,
                            command=self.tree.yview)
        vsb.pack(side=tk.RIGHT, fill=tk.Y)
        self.tree.configure(yscrollcommand=vsb.set)
        self.tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        self.tree.bind("<Double-Button-1>", self._on_double_click)

        # Bottom buttons
        btn = ttk.Frame(self.window, padding=6)
        btn.pack(fill=tk.X)
        ttk.Button(btn, text="Apply",  command=self._apply_to_globals)\
            .pack(side=tk.LEFT, padx=4)
        ttk.Button(btn, text="Reload", command=self.load_data)\
            .pack(side=tk.LEFT, padx=4)
        ttk.Button(btn, text="Close",  command=self._on_close)\
            .pack(side=tk.RIGHT, padx=4)

    def load_data(self):
        """Load current AWG waveguide data from globals."""
        if not self.window or not self.window.winfo_exists():
            return
        self.tree.delete(*self.tree.get_children())
        for j in range(1, g.na + 1):
            self.tree.insert("", "end", iid=str(j), values=(
                j,
                f"{g.Awg_Loss[j]:.2f}",
                f"{g.Awg_Phase[j]:.2f}",
                f"{g.Awg_Length[j]:.2f}",
                f"{g.Awg_Arc_R[j]:.2f}",
                f"{g.Awg_Arc_L[j]:.2f}",
                f"{g.Awg_Width[j]:.2f}",
                f"{g.Awg_Height[j]:.2f}",
            ))

    def _apply_to_globals(self):
        """Write edited values back to globals."""
        for item in self.tree.get_children():
            vals = self.tree.item(item, "values")
            try:
                j = int(vals[0])
                if 1 <= j <= g.na:
                    g.Awg_Loss[j]   = float(vals[1])
                    g.Awg_Phase[j]  = float(vals[2])
                    g.Awg_Length[j] = float(vals[3])
                    g.Awg_Arc_R[j]  = float(vals[4])
                    g.Awg_Arc_L[j]  = float(vals[5])
                    g.Awg_Width[j]  = float(vals[6])
                    g.Awg_Height[j] = float(vals[7])
            except (ValueError, IndexError):
                pass

    def _apply_loss(self):
        try:
            loss = float(self.var_loss.get())
        except ValueError:
            return
        for item in self.tree.get_children():
            vals = list(self.tree.item(item, "values"))
            vals[1] = f"{loss:.2f}"
            self.tree.item(item, values=vals)
        for j in range(1, g.na + 1):
            g.Awg_Loss[j] = loss

    def _zero_phase(self):
        for item in self.tree.get_children():
            vals = list(self.tree.item(item, "values"))
            vals[2] = "0.00"
            self.tree.item(item, values=vals)
        for j in range(1, g.na + 1):
            g.Awg_Phase[j] = 0.0

    def _random_phase(self):
        import random
        for item in self.tree.get_children():
            vals = list(self.tree.item(item, "values"))
            p = random.uniform(-180, 180)
            vals[2] = f"{p:.2f}"
            self.tree.item(item, values=vals)
        self._apply_to_globals()

    def _on_double_click(self, event):
        """Inline edit cell on double-click."""
        region = self.tree.identify_region(event.x, event.y)
        if region != "cell":
            return
        item = self.tree.identify_row(event.y)
        col  = self.tree.identify_column(event.x)

        # Map column number to index (col is "#1", "#2", ...)
        col_idx = int(col.replace("#", "")) - 1
        if col_idx == 0:  # row number, not editable
            return

        x, y, w, h = self.tree.bbox(item, col)
        vals = list(self.tree.item(item, "values"))
        entry_var = tk.StringVar(value=vals[col_idx])

        entry = ttk.Entry(self.tree, textvariable=entry_var, width=10)
        entry.place(x=x, y=y, width=w, height=h)
        entry.focus()

        def _confirm(event=None):
            vals[col_idx] = entry_var.get()
            self.tree.item(item, values=vals)
            entry.destroy()

        entry.bind("<Return>",    _confirm)
        entry.bind("<FocusOut>",  _confirm)
        entry.bind("<Escape>",    lambda e: entry.destroy())

    def _on_close(self):
        self._apply_to_globals()
        if self.window:
            self.window.withdraw()
