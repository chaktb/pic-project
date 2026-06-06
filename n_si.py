"""Refractive index of crystalline silicon (Si) via the Sellmeier equation.

Salzberg & Villa (1957), as tabulated for crystalline silicon:

    n(lambda)^2 = 1 + B1*l^2/(l^2 - C1) + B2*l^2/(l^2 - C2) + B3*l^2/(l^2 - C3)

    B1 = 10.6684293,    C1 = 0.301516485^2
    B2 = 0.0030434748,  C2 = 1.13475115^2
    B3 = 1.54133408,    C3 = 1104.0^2

Valid range: 1.36 - 11 um.

Usage:
    python n_si.py                 # save n_si.png (index vs wavelength)
    python n_si.py 1.55            # print n and n_g at 1.55 um
    python n_si.py 1.31 1.55 3.0   # several wavelengths
"""

import math
import sys

B1, B2, B3 = 10.6684293, 0.0030434748, 1.54133408
C1, C2, C3 = 0.301516485**2, 1.13475115**2, 1104.0**2


def n_si(wavelength_um):
    """Refractive index of crystalline silicon via Sellmeier (Salzberg & Villa).

    Valid range: 1.36 - 11 um. Takes a wavelength in micrometres (scalar).
    """
    l2 = wavelength_um**2
    n2 = 1 + (B1 * l2) / (l2 - C1) + (B2 * l2) / (l2 - C2) + (B3 * l2) / (l2 - C3)
    return math.sqrt(n2)


def n_group(wavelength_um):
    """Group index of silicon, n_g = n - lambda * dn/dlambda.

    Differentiating the Sellmeier equation gives the closed form
        n_g = n + (lambda^2 / n) * sum_i B_i * C_i / (lambda^2 - C_i)^2.
    """
    l2 = wavelength_um**2
    n = n_si(wavelength_um)
    s = (B1 * C1 / (l2 - C1) ** 2
         + B2 * C2 / (l2 - C2) ** 2
         + B3 * C3 / (l2 - C3) ** 2)
    return n + (l2 / n) * s


def _save_plot():
    import numpy as np
    import matplotlib.pyplot as plt

    wavelengths = np.linspace(1.4, 6.0, 1000)
    n = np.array([n_si(w) for w in wavelengths])
    ng = np.array([n_group(w) for w in wavelengths])

    plt.figure(figsize=(8, 5))
    plt.plot(wavelengths, n, color="tab:blue", linewidth=2, label="n (phase index)")
    plt.plot(wavelengths, ng, color="tab:orange", linewidth=2, label="n_g (group index)")
    plt.xlabel("Wavelength (μm)")
    plt.ylabel("Index")
    plt.title("Silicon: Phase and Group Index (Sellmeier, Salzberg & Villa)")
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.legend()
    plt.tight_layout()
    plt.savefig("n_si.png", dpi=150)
    print("Saved: n_si.png")


if __name__ == "__main__":
    if len(sys.argv) == 1:
        _save_plot()
    else:
        for arg in sys.argv[1:]:
            wl = float(arg)
            print(f"lambda = {wl:.4f} um   n = {n_si(wl):.6f}   n_g = {n_group(wl):.6f}")
