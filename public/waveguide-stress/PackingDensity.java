/**
 * PECVD Ge-doped SiO2 packing density — reference implementation.
 *
 * Three estimators (matching the web page):
 *   A. Lorentz-Lorentz    : from measured refractive index n
 *   B. Density ratio       : p = rhoFilm / rhoBulk
 *   C. Deposition-condition: empirical Arrhenius/oxidation/ion model
 */
public class PackingDensity {

    /** Fused silica Sellmeier index (Malitson 1965). */
    static double nSiO2(double wlUm) {
        double l2 = wlUm * wlUm;
        double n2 = 1
                + 0.6961663 * l2 / (l2 - Math.pow(0.0684043, 2))
                + 0.4079426 * l2 / (l2 - Math.pow(0.1162414, 2))
                + 0.8974794 * l2 / (l2 - Math.pow(9.896161, 2));
        return Math.sqrt(n2);
    }

    /** GeO2 glass Sellmeier index (Fleming 1984). */
    static double nGeO2(double wlUm) {
        double l2 = wlUm * wlUm;
        double n2 = 1
                + 0.80686642 * l2 / (l2 - Math.pow(0.068972606, 2))
                + 0.71815848 * l2 / (l2 - Math.pow(0.15396605, 2))
                + 0.85416831 * l2 / (l2 - Math.pow(11.841931, 2));
        return Math.sqrt(n2);
    }

    /** Lorentz-Lorentz molar refraction term. */
    static double lorentz(double n) {
        return (n * n - 1) / (n * n + 2);
    }

    /** Molar-averaged bulk index of (1-x) SiO2 + x GeO2. */
    static double nBulk(double xGeO2, double wlUm) {
        double Rm = (1 - xGeO2) * lorentz(nSiO2(wlUm)) + xGeO2 * lorentz(nGeO2(wlUm));
        return Math.sqrt((1 + 2 * Rm) / (1 - Rm));
    }

    /** Bulk density [g/cm3], linear in mol fraction. */
    static double rhoBulk(double xGeO2) {
        return (1 - xGeO2) * 2.20 + xGeO2 * 3.65;
    }

    // ---- Method A ----
    static double packingFromIndex(double nFilm, double xGeO2, double wlUm) {
        return lorentz(nFilm) / lorentz(nBulk(xGeO2, wlUm));
    }

    // ---- Method B ----
    static double packingFromDensity(double rhoFilm, double xGeO2) {
        return rhoFilm / rhoBulk(xGeO2);
    }

    // ---- Method C ----
    static double packingFromDeposition(double sih4, double n2o, double geh4,
                                        double tempC, double pTorr, double rfW) {
        double p0 = 0.94, Ea = 0.15, k = 8.617e-5, Tref = 573.15;
        double Tk = tempC + 273.15;
        double fT = Math.exp(-(Ea / k) * (1 / Tk - 1 / Tref));
        double rOx = n2o / sih4;
        double fOx = 1 - 0.08 * Math.exp(-rOx / 15);
        double fIon = 1 + 0.015 * Math.log(1 + (rfW / pTorr) / 150);
        return Math.min(p0 * fT * fOx * fIon, 1.02);
    }

    public static void main(String[] args) {
        double pA = packingFromIndex(1.448, 0.05, 1.55);
        System.out.printf("A. Lorentz-Lorentz : n=1.448 x=0.05  n_bulk=%.4f  p=%.4f%n",
                nBulk(0.05, 1.55), pA);

        double pB = packingFromDensity(2.15, 0.05);
        System.out.printf("B. Density ratio    : rho=2.15 x=0.05  rho_bulk=%.3f  p=%.4f%n",
                rhoBulk(0.05), pB);

        double pC = packingFromDeposition(20, 900, 5, 300, 1.0, 150);
        System.out.printf("C. Deposition        : T=300C  p=%.4f%n", pC);
    }
}
