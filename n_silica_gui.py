"""Tkinter GUI: enter a wavelength (um) and read the fused-silica refractive index.

Uses the Sellmeier equation (Malitson 1965) from n_silica.py.

Run: python n_silica_gui.py
"""

import tkinter as tk
from tkinter import ttk

from n_silica import n_silica

VALID_MIN, VALID_MAX = 0.21, 6.7


def main():
    root = tk.Tk()
    root.title("Refractive Index of Fused Silica")
    root.resizable(False, False)

    frame = ttk.Frame(root, padding=16)
    frame.grid(sticky="nsew")

    ttk.Label(frame, text="Sellmeier equation (Malitson 1965)",
              font=("TkDefaultFont", 10, "bold")).grid(
        row=0, column=0, columnspan=3, sticky="w", pady=(0, 12))

    ttk.Label(frame, text="Wavelength λ").grid(row=1, column=0, sticky="w")
    wl_var = tk.StringVar(value="1.55")
    entry = ttk.Entry(frame, textvariable=wl_var, width=12, justify="right")
    entry.grid(row=1, column=1, padx=6)
    ttk.Label(frame, text="μm").grid(row=1, column=2, sticky="w")

    ttk.Separator(frame, orient="horizontal").grid(
        row=2, column=0, columnspan=3, sticky="ew", pady=12)

    ttk.Label(frame, text="Refractive index n").grid(row=3, column=0, sticky="w")
    result_var = tk.StringVar(value="—")
    ttk.Label(frame, textvariable=result_var,
              font=("TkDefaultFont", 14, "bold")).grid(
        row=3, column=1, columnspan=2, sticky="e")

    status_var = tk.StringVar(value=f"Valid range: {VALID_MIN} – {VALID_MAX} μm")
    status = ttk.Label(frame, textvariable=status_var, foreground="#888")
    status.grid(row=4, column=0, columnspan=3, sticky="w", pady=(12, 0))

    def compute(*_):
        text = wl_var.get().strip()
        try:
            wl = float(text)
        except ValueError:
            result_var.set("—")
            status_var.set("Enter a number." if text else
                           f"Valid range: {VALID_MIN} – {VALID_MAX} μm")
            return
        result_var.set(f"{float(n_silica(wl)):.6f}")
        if wl < VALID_MIN or wl > VALID_MAX:
            status_var.set(f"⚠ Outside valid range ({VALID_MIN} – {VALID_MAX} μm)")
        else:
            status_var.set(f"Valid range: {VALID_MIN} – {VALID_MAX} μm")

    wl_var.trace_add("write", compute)
    entry.focus_set()
    entry.selection_range(0, "end")
    compute()

    root.mainloop()


if __name__ == "__main__":
    main()
