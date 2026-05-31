"""
About Dialog
Ported from Frm_ID.frm (VB 5.0/6.0) to Python/Tkinter
"""

import tkinter as tk
from tkinter import ttk
import globals as g


class FrmAbout:
    """About dialog box."""

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
        self.window.title("About")
        self.window.resizable(False, False)
        self.window.protocol("WM_DELETE_WINDOW", self._on_ok)

        # Center it
        self.window.geometry("+%d+%d" % (
            self.master.winfo_screenwidth() // 2 - 200,
            self.master.winfo_screenheight() // 2 - 100))

        ttk.Label(self.window,
                  text=g.VER_STRING,
                  font=("Arial", 13, "bold"))\
            .pack(padx=24, pady=(20, 6))

        ttk.Label(self.window,
                  text=g.COPYRIGHT,
                  font=("Arial", 10))\
            .pack(padx=24, pady=4)

        ttk.Label(self.window,
                  text="Python/Tkinter port – 2024",
                  font=("Arial", 9), foreground="#666666")\
            .pack(padx=24, pady=(4, 16))

        ttk.Button(self.window, text="OK", width=10, command=self._on_ok)\
            .pack(pady=(0, 16))

    def _on_ok(self):
        if self.window:
            self.window.withdraw()
