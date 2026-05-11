import numpy as np
import matplotlib.pyplot as plt


def n_silica(wavelength_um):
    """Refractive index of fused silica via Sellmeier equation (Malitson 1965).

    Valid range: 0.21 - 6.7 um.
    """
    B1, B2, B3 = 0.6961663, 0.4079426, 0.8974794
    C1, C2, C3 = 0.0684043**2, 0.1162414**2, 9.896161**2

    l2 = wavelength_um**2
    n2 = 1 + (B1 * l2) / (l2 - C1) + (B2 * l2) / (l2 - C2) + (B3 * l2) / (l2 - C3)
    return np.sqrt(n2)


if __name__ == "__main__":
    wavelengths = np.linspace(0.25, 2.5, 1000)
    n = n_silica(wavelengths)

    plt.figure(figsize=(8, 5))
    plt.plot(wavelengths, n, color="tab:blue", linewidth=2)
    plt.xlabel("Wavelength (μm)")
    plt.ylabel("Refractive index n")
    plt.title("Refractive Index of Fused Silica (Sellmeier Equation)")
    plt.grid(True, linestyle="--", alpha=0.6)
#    plt.tight_layout()
#    plt.show()
    
    plt.tight_layout()
    plt.savefig("n_silica.png", dpi=150)
    print("Saved: n_silica.png")
    