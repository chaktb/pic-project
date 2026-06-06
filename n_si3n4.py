"""Refractive index of silicon nitride (Si3N4) via the Sellmeier equation.

Luke et al., Opt. Lett. 40, 4823 (2015) — LPCVD stoichiometric Si3N4:

    n(lambda)^2 = 1 + B1*lambda^2/(lambda^2 - C1) + B2*lambda^2/(lambda^2 - C2)

    B1 = 3.0249,   C1 = 0.1353406^2
    B2 = 40314.0,  C2 = 1239.842^2

Valid range: 0.31 - 5.504 um.

Usage:
    python n_si3n4.py                 # save n_si3n4.png (index vs wavelength)
    python n_si3n4.py 1.55            # print n and n_g at 1.55 um
    python n_si3n4.py 1.31 1.55 0.85  # several wavelengths
"""

import math
import sys

B1, B2 = 3.0249, 40314.0
C1, C2 = 0.1353406**2, 1239.842**2


def n_si3n4(wavelength_um):
    """Refractive index of Si3N4 via Sellmeier equation (Luke 2015).

    Valid range: 0.31 - 5.504 um. Takes a wavelength in micrometres (scalar).
    """
    l2 = wavelength_um**2
    n2 = 1 + (B1 * l2) / (l2 - C1) + (B2 * l2) / (l2 - C2)
    return math.sqrt(n2)


def n_group(wavelength_um):
    """Group index of Si3N4, n_g = n - lambda * dn/dlambda.

    Differentiating the Sellmeier equation gives the closed form
        n_g = n + (lambda^2 / n) * sum_i B_i * C_i / (lambda^2 - C_i)^2.
    """
    l2 = wavelength_um**2
    n = n_si3n4(wavelength_um)
    s = B1 * C1 / (l2 - C1) ** 2 + B2 * C2 / (l2 - C2) ** 2
    return n + (l2 / n) * s


def _save_plot():
    import numpy as np
    import matplotlib.pyplot as plt

    wavelengths = np.linspace(0.4, 2.5, 1000)
    n = np.array([n_si3n4(w) for w in wavelengths])
    ng = np.array([n_group(w) for w in wavelengths])

    plt.figure(figsize=(8, 5))
    plt.plot(wavelengths, n, color="tab:blue", linewidth=2, label="n (phase index)")
    plt.plot(wavelengths, ng, color="tab:orange", linewidth=2, label="n_g (group index)")
    plt.xlabel("Wavelength (μm)")
    plt.ylabel("Index")
    plt.title("Si₃N₄: Phase and Group Index (Sellmeier, Luke 2015)")
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.legend()
    plt.tight_layout()
    plt.savefig("n_si3n4.png", dpi=150)
    print("Saved: n_si3n4.png")


if __name__ == "__main__":
    if len(sys.argv) == 1:
        _save_plot()
    else:
        for arg in sys.argv[1:]:
            wl = float(arg)
            print(f"lambda = {wl:.4f} um   n = {n_si3n4(wl):.6f}   n_g = {n_group(wl):.6f}")
