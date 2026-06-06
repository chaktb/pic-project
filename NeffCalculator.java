/**
 * Effective index of a silica strip waveguide via EIM (TE-TE).
 *
 * Inputs: wavelength (um), core width W (um), core height H (um),
 *         index contrast delta = (n_core - n_clad) / n_core * 100 (%).
 *
 * n_clad is the fused-silica Sellmeier index (Malitson 1965);
 * n_core is derived from the contrast: n_core = n_clad / (1 - delta/100).
 *
 * Compile: javac NeffCalculator.java
 * Run:     java NeffCalculator 1.55 5.0 5.0 0.75
 */
public class NeffCalculator {

    public static double nSilica(double wavelengthUm) {
        double B1 = 0.6961663, B2 = 0.4079426, B3 = 0.8974794;
        double C1 = 0.0684043 * 0.0684043;
        double C2 = 0.1162414 * 0.1162414;
        double C3 = 9.896161 * 9.896161;
        double l2 = wavelengthUm * wavelengthUm;
        double n2 = 1.0 + B1 * l2 / (l2 - C1) + B2 * l2 / (l2 - C2) + B3 * l2 / (l2 - C3);
        return Math.sqrt(n2);
    }

    public static double slabNeffTE(double wlUm, double dUm, double nCore, double nClad) {
        if (nCore <= nClad) return nClad;
        double k0 = 2.0 * Math.PI / wlUm;
        double V = (k0 * dUm / 2.0) * Math.sqrt(nCore * nCore - nClad * nClad);
        if (V == 0.0) return nClad;

        double upper = Math.min(Math.PI / 2.0 - 1e-9, V - 1e-12);
        if (upper <= 0.0) return nClad;

        double lo = 1e-9, hi = upper;
        double flo = f(lo, V), fhi = f(hi, V);
        if (flo * fhi > 0.0) return nClad;

        for (int i = 0; i < 100; i++) {
            double mid = 0.5 * (lo + hi);
            double fmid = f(mid, V);
            if (flo * fmid <= 0.0) {
                hi = mid; fhi = fmid;
            } else {
                lo = mid; flo = fmid;
            }
            if (hi - lo < 1e-12) break;
        }
        double u = 0.5 * (lo + hi);
        double kx = 2.0 * u / dUm;
        double beta2 = (nCore * k0) * (nCore * k0) - kx * kx;
        return Math.sqrt(beta2) / k0;
    }

    private static double f(double u, double V) {
        return u * Math.tan(u) - Math.sqrt(Math.max(V * V - u * u, 0.0));
    }

    public static double[] computeNeff(double wlUm, double wUm, double hUm, double deltaPct) {
        double nClad = nSilica(wlUm);
        double nCore = nClad / (1.0 - deltaPct / 100.0);
        double nEffV = slabNeffTE(wlUm, hUm, nCore, nClad);
        double nEff  = slabNeffTE(wlUm, wUm, nEffV, nClad);
        return new double[] { nCore, nClad, nEffV, nEff };
    }

    public static void main(String[] args) {
        if (args.length != 4) {
            System.out.println("Usage: java NeffCalculator <wavelength_um> <width_um> <height_um> <delta_percent>");
            System.out.println("Example: java NeffCalculator 1.55 5.0 5.0 0.75");
            return;
        }
        double wl = Double.parseDouble(args[0]);
        double w  = Double.parseDouble(args[1]);
        double h  = Double.parseDouble(args[2]);
        double d  = Double.parseDouble(args[3]);
        double[] r = computeNeff(wl, w, h, d);
        System.out.printf("wavelength            = %.4f um%n", wl);
        System.out.printf("core size  W x H      = %.4f x %.4f um%n", w, h);
        System.out.printf("index contrast        = %.4f %%%n", d);
        System.out.printf("n_clad (Sellmeier)    = %.6f%n", r[1]);
        System.out.printf("n_core (from delta)   = %.6f%n", r[0]);
        System.out.printf("n_eff (vertical slab) = %.6f%n", r[2]);
        System.out.printf("n_eff (full EIM)      = %.6f%n", r[3]);
    }
}
