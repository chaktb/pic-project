"""Refractive index of fused silica via the Sellmeier equation (Malitson 1965).

n(lambda)^2 = 1 + sum_i  B_i * lambda^2 / (lambda^2 - C_i)

    B1 = 0.6961663,  C1 = 0.0684043^2
    B2 = 0.4079426,  C2 = 0.1162414^2
    B3 = 0.8974794,  C3 = 9.896161^2

Valid range: 0.21 - 6.7 um.

Usage:
    python n_silica.py                 # save n_silica.png (index vs wavelength)
    python n_silica.py 1.55            # print n at 1.55 um
    python n_silica.py 1.31 1.55 0.85  # print n at several wavelengths
"""

import math
import sys


def n_silica(wavelength_um):
    """Refractive index of fused silica via Sellmeier equation (Malitson 1965).

    Valid range: 0.21 - 6.7 um. Takes a wavelength in micrometres (scalar).
    """
    B1, B2, B3 = 0.6961663, 0.4079426, 0.8974794
    C1, C2, C3 = 0.0684043**2, 0.1162414**2, 9.896161**2

    l2 = wavelength_um**2
    n2 = 1 + (B1 * l2) / (l2 - C1) + (B2 * l2) / (l2 - C2) + (B3 * l2) / (l2 - C3)
    return math.sqrt(n2)


def _save_plot():
    import numpy as np
    import matplotlib.pyplot as plt

    wavelengths = np.linspace(0.25, 2.5, 1000)
    n = np.array([n_silica(w) for w in wavelengths])

    plt.figure(figsize=(8, 5))
    plt.plot(wavelengths, n, color="tab:blue", linewidth=2)
    plt.xlabel("Wavelength (μm)")
    plt.ylabel("Refractive index n")
    plt.title("Refractive Index of Fused Silica (Sellmeier Equation)")
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.tight_layout()
    plt.savefig("n_silica.png", dpi=150)
    print("Saved: n_silica.png")


if __name__ == "__main__":
    if len(sys.argv) == 1:
        _save_plot()
    else:
        for arg in sys.argv[1:]:
            wl = float(arg)
            print(f"lambda = {wl:.4f} um   n = {n_silica(wl):.6f}")
