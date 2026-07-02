#!/usr/bin/env python3
"""PECVD Ge-doped SiO2 packing density — reference implementation.

Three independent estimators (matching the web page):
  A. Lorentz-Lorentz    : from measured refractive index n
  B. Density ratio       : p = rho_film / rho_bulk
  C. Deposition-condition: empirical Arrhenius/oxidation/ion model

Packing density p = (film molar density) / (fully-dense bulk).
p < 1 indicates porosity / voids / incomplete network (e.g. Si-OH, Si-H).
"""

import math


# ---------- Sellmeier bulk indices ----------
def n_sio2(wl_um: float) -> float:
    """Fused silica, Malitson 1965."""
    l2 = wl_um * wl_um
    n2 = (1
          + 0.6961663 * l2 / (l2 - 0.0684043**2)
          + 0.4079426 * l2 / (l2 - 0.1162414**2)
          + 0.8974794 * l2 / (l2 - 9.896161**2))
    return math.sqrt(n2)


def n_geo2(wl_um: float) -> float:
    """GeO2 glass, Fleming 1984."""
    l2 = wl_um * wl_um
    n2 = (1
          + 0.80686642 * l2 / (l2 - 0.068972606**2)
          + 0.71815848 * l2 / (l2 - 0.15396605**2)
          + 0.85416831 * l2 / (l2 - 11.841931**2))
    return math.sqrt(n2)


def _R(n: float) -> float:
    """Lorentz-Lorentz molar refraction term."""
    return (n * n - 1) / (n * n + 2)


def n_bulk(x_geo2: float, wl_um: float) -> float:
    """Molar-averaged bulk index of (1-x) SiO2 + x GeO2."""
    Rm = (1 - x_geo2) * _R(n_sio2(wl_um)) + x_geo2 * _R(n_geo2(wl_um))
    return math.sqrt((1 + 2 * Rm) / (1 - Rm))


def rho_bulk(x_geo2: float) -> float:
    """Bulk density [g/cm3], linear in mol fraction (SiO2 2.20, GeO2 3.65)."""
    return (1 - x_geo2) * 2.20 + x_geo2 * 3.65


# ---------- Method A: Lorentz-Lorentz ----------
def packing_from_index(n_film: float, x_geo2: float = 0.0, wl_um: float = 1.55) -> float:
    nb = n_bulk(x_geo2, wl_um)
    return _R(n_film) / _R(nb)


# ---------- Method B: density ratio ----------
def packing_from_density(rho_film: float, x_geo2: float = 0.0) -> float:
    return rho_film / rho_bulk(x_geo2)


# ---------- Method C: deposition-condition empirical ----------
def packing_from_deposition(sih4_sccm: float, n2o_sccm: float, geh4_sccm: float,
                            temp_c: float, pressure_torr: float, rf_w: float):
    """Empirical densification model. CALIBRATE constants to your reactor.

    p = p0 * f_T * f_ox * f_ion
      f_T   = exp[-(Ea/k)(1/T - 1/Tref)]     temperature (Arrhenius)
      f_ox  = 1 - A*exp(-r_ox/r0)            oxidant ratio N2O/SiH4
      f_ion = 1 + B*ln(1 + (RF/P)/S0)        ion-bombardment energy
    """
    p0, Ea, k, Tref = 0.94, 0.15, 8.617e-5, 573.15  # 300 C reference
    Tk = temp_c + 273.15
    f_T = math.exp(-(Ea / k) * (1 / Tk - 1 / Tref))
    r_ox = n2o_sccm / sih4_sccm
    f_ox = 1 - 0.08 * math.exp(-r_ox / 15)
    f_ion = 1 + 0.015 * math.log(1 + (rf_w / pressure_torr) / 150)
    p = min(p0 * f_T * f_ox * f_ion, 1.02)
    x_ge = 0.6 * geh4_sccm / (sih4_sccm + geh4_sccm)   # incorporation eff ~0.6
    return {"p": p, "f_T": f_T, "f_ox": f_ox, "f_ion": f_ion,
            "r_ox": r_ox, "x_ge": x_ge}


if __name__ == "__main__":
    print("A. Lorentz-Lorentz : n=1.470, x=0.05, 1.55um")
    pA = packing_from_index(1.470, 0.05, 1.55)
    print(f"   n_bulk={n_bulk(0.05,1.55):.4f}  p={pA:.4f}  void={(1-pA)*100:.2f}%\n")

    print("B. Density ratio    : rho_film=2.15, x=0.05")
    pB = packing_from_density(2.15, 0.05)
    print(f"   rho_bulk={rho_bulk(0.05):.3f}  p={pB:.4f}  void={(1-pB)*100:.2f}%\n")

    print("C. Deposition        : SiH4=20 N2O=900 GeH4=5 T=300C P=1.0Torr RF=150W")
    r = packing_from_deposition(20, 900, 5, 300, 1.0, 150)
    print(f"   r_ox={r['r_ox']:.1f}  f_T={r['f_T']:.4f}  f_ox={r['f_ox']:.4f}  "
          f"f_ion={r['f_ion']:.4f}")
    print(f"   x_Ge={r['x_ge']:.4f}  p={r['p']:.4f}  void={(1-r['p'])*100:.2f}%")
