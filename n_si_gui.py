"""Tkinter GUI: enter a wavelength (um) and read the silicon phase index n
and group index n_g, with both points marked on a live Sellmeier plot.

Uses the Sellmeier equation (Salzberg & Villa) from n_si.py.

Requires: matplotlib (for the embedded graph).
Run: python n_si_gui.py
"""

import tkinter as tk
from tkinter import ttk

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from n_si import n_si, n_group

VALID_MIN, VALID_MAX = 1.36, 11.0  # Salzberg & Villa validity
PLOT_MIN, PLOT_MAX = 1.4, 6.0      # plotted range


def main():
    root = tk.Tk()
    root.title("Refractive Index of Si")

    frame = ttk.Frame(root, padding=12)
    frame.grid(sticky="nsew")

    # --- input row ---
    ttk.Label(frame, text="Wavelength λ").grid(row=0, column=0, sticky="w")
    wl_var = tk.StringVar(value="1.55")
    entry = ttk.Entry(frame, textvariable=wl_var, width=12, justify="right")
    entry.grid(row=0, column=1, padx=6)
    ttk.Label(frame, text="μm").grid(row=0, column=2, sticky="w")

    ttk.Label(frame, text="n =").grid(row=0, column=3, padx=(18, 4), sticky="e")
    n_var = tk.StringVar(value="—")
    ttk.Label(frame, textvariable=n_var, foreground="#1f77b4",
              font=("TkDefaultFont", 14, "bold")).grid(row=0, column=4, sticky="w")
    ttk.Label(frame, text="n_g =").grid(row=0, column=5, padx=(18, 4), sticky="e")
    ng_var = tk.StringVar(value="—")
    ttk.Label(frame, textvariable=ng_var, foreground="#e8590c",
              font=("TkDefaultFont", 14, "bold")).grid(row=0, column=6, sticky="w")

    status_var = tk.StringVar(value=f"Valid range: {VALID_MIN} – {VALID_MAX} μm")
    ttk.Label(frame, textvariable=status_var, foreground="#888").grid(
        row=1, column=0, columnspan=7, sticky="w", pady=(6, 10))

    # --- plot ---
    fig = Figure(figsize=(6.8, 4.0), dpi=100)
    ax = fig.add_subplot(111)

    wls = [PLOT_MIN + (PLOT_MAX - PLOT_MIN) * i / 400 for i in range(401)]
    ax.plot(wls, [n_si(w) for w in wls], color="tab:blue", linewidth=2, label="n (phase index)")
    ax.plot(wls, [n_group(w) for w in wls], color="tab:orange", linewidth=2, label="n_g (group index)")
    ax.set_xlabel("Wavelength λ (μm)")
    ax.set_ylabel("Index")
    ax.set_title("Silicon: Phase and Group Index (Sellmeier, Salzberg & Villa)")
    ax.grid(True, linestyle="--", alpha=0.6)
    ax.set_xlim(PLOT_MIN, PLOT_MAX)
    ax.legend(loc="upper right")

    marker_n, = ax.plot([], [], "o", color="tab:blue", markersize=8, zorder=5)
    marker_g, = ax.plot([], [], "o", color="tab:orange", markersize=8, zorder=5)
    guide_v = ax.axvline(x=0, color="gray", linestyle=":", linewidth=1, visible=False)
    annot_n = ax.annotate("", xy=(0, 0), xytext=(8, -14), textcoords="offset points",
                          fontsize=9, color="white",
                          bbox=dict(boxstyle="round,pad=0.3", fc="tab:blue", ec="none"),
                          visible=False)
    annot_g = ax.annotate("", xy=(0, 0), xytext=(8, 8), textcoords="offset points",
                          fontsize=9, color="white",
                          bbox=dict(boxstyle="round,pad=0.3", fc="tab:orange", ec="none"),
                          visible=False)
    fig.tight_layout()

    canvas = FigureCanvasTkAgg(fig, master=frame)
    canvas.get_tk_widget().grid(row=2, column=0, columnspan=7, sticky="nsew")

    artists = [marker_n, marker_g, guide_v, annot_n, annot_g]

    def set_visible(v):
        for a in artists:
            a.set_visible(v)

    def compute(*_):
        text = wl_var.get().strip()
        try:
            wl = float(text)
        except ValueError:
            n_var.set("—"); ng_var.set("—")
            status_var.set("Enter a number." if text else
                           f"Valid range: {VALID_MIN} – {VALID_MAX} μm")
            set_visible(False)
            canvas.draw_idle()
            return

        n = n_si(wl)
        g = n_group(wl)
        n_var.set(f"{n:.6f}")
        ng_var.set(f"{g:.6f}")

        if wl < VALID_MIN or wl > VALID_MAX:
            status_var.set(f"⚠ Outside Sellmeier validity ({VALID_MIN} – {VALID_MAX} μm)")
        elif wl < PLOT_MIN or wl > PLOT_MAX:
            status_var.set(f"values shown; markers pinned to plotted range {PLOT_MIN} – {PLOT_MAX} μm")
        else:
            status_var.set(f"λ = {wl} μm   →   n = {n:.6f},   n_g = {g:.6f}")

        wl_c = min(PLOT_MAX, max(PLOT_MIN, wl))
        n_c, g_c = n_si(wl_c), n_group(wl_c)
        marker_n.set_data([wl_c], [n_c])
        marker_g.set_data([wl_c], [g_c])
        guide_v.set_xdata([wl_c, wl_c])
        annot_n.xy = (wl_c, n_c); annot_n.set_text(f"n={n_c:.4f}")
        annot_g.xy = (wl_c, g_c); annot_g.set_text(f"n_g={g_c:.4f}")
        set_visible(True)
        canvas.draw_idle()

    wl_var.trace_add("write", compute)
    entry.focus_set()
    entry.selection_range(0, "end")
    compute()

    root.mainloop()


if __name__ == "__main__":
    main()
