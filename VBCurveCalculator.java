/**
 * V-b curve table for a symmetric step-index slab waveguide.
 *
 * V = k0 * a * sqrt(n_core^2 - n_clad^2),  a = d/2
 * b = (n_eff^2 - n_clad^2) / (n_core^2 - n_clad^2)
 * TE_m cutoff at V_c = m * pi/2.  Single mode: V < pi/2.
 *
 * Compile: javac VBCurveCalculator.java
 * Run:     java VBCurveCalculator 1.55 5.0 0.75
 */
public class VBCurveCalculator {

    public static double nSilica(double wlUm) {
        double B1 = 0.6961663, B2 = 0.4079426, B3 = 0.8974794;
        double C1 = 0.0684043 * 0.0684043;
        double C2 = 0.1162414 * 0.1162414;
        double C3 = 9.896161 * 9.896161;
        double l2 = wlUm * wlUm;
        double n2 = 1.0 + B1 * l2 / (l2 - C1) + B2 * l2 / (l2 - C2) + B3 * l2 / (l2 - C3);
        return Math.sqrt(n2);
    }

    /** b for TE_m at given V; returns NaN if below cutoff. */
    public static double bForMode(double V, int m) {
        double Vc = m * Math.PI / 2.0;
        if (V <= Vc) return Double.NaN;
        boolean even = (m % 2 == 0);
        double uMin = Math.max(Vc + 1e-9, 1e-9);
        double uMax = Math.min((m + 1) * Math.PI / 2.0 - 1e-9, V - 1e-12);
        if (uMax <= uMin) return Double.NaN;

        double lo = uMin, hi = uMax;
        double flo = f(lo, V, even), fhi = f(hi, V, even);
        if (flo * fhi > 0) return Double.NaN;

        for (int i = 0; i < 100; i++) {
            double mid = 0.5 * (lo + hi);
            double fmid = f(mid, V, even);
            if (flo * fmid <= 0) { hi = mid; fhi = fmid; }
            else { lo = mid; flo = fmid; }
            if (hi - lo < 1e-12) break;
        }
        double u = 0.5 * (lo + hi);
        double w = Math.sqrt(Math.max(V * V - u * u, 0.0));
        return (w / V) * (w / V);
    }

    private static double f(double u, double V, boolean even) {
        double w = Math.sqrt(Math.max(V * V - u * u, 0.0));
        return even ? (u * Math.tan(u) - w) : (-u / Math.tan(u) - w);
    }

    public static double operatingV(double wlUm, double dUm, double nCore, double nClad) {
        if (nCore <= nClad) return 0.0;
        double k0 = 2.0 * Math.PI / wlUm;
        double a = dUm / 2.0;
        return k0 * a * Math.sqrt(nCore * nCore - nClad * nClad);
    }

    public static void main(String[] args) {
        if (args.length != 3) {
            System.out.println("Usage: java VBCurveCalculator <wavelength_um> <thickness_um> <delta_percent>");
            return;
        }
        double wl = Double.parseDouble(args[0]);
        double d  = Double.parseDouble(args[1]);
        double dl = Double.parseDouble(args[2]);

        double nCore = nSilica(wl);
        double nClad = nCore * (1.0 - dl / 100.0);
        double Vop = operatingV(wl, d, nCore, nClad);
        boolean single = Vop < Math.PI / 2.0;

        System.out.printf("wavelength      = %.4f um%n", wl);
        System.out.printf("core thickness  = %.4f um%n", d);
        System.out.printf("index contrast  = %.4f %%%n", dl);
        System.out.printf("n_core (Sellmeier) = %.6f%n", nCore);
        System.out.printf("n_clad             = %.6f%n", nClad);
        System.out.printf("V_operating        = %.4f%n", Vop);
        System.out.printf("V_cutoff (TE_1)    = %.4f%n", Math.PI / 2.0);
        System.out.printf("single mode        = %s%n", single ? "YES" : "NO");
        System.out.println();

        System.out.println("V-b table (b = NaN means mode below cutoff):");
        System.out.println("    V        b(TE_0)  b(TE_1)  b(TE_2)  b(TE_3)");
        double Vmax = Math.max(10.0, Vop + 2);
        for (double V = 0.5; V <= Vmax + 1e-9; V += 0.5) {
            System.out.printf("  %5.2f", V);
            for (int m = 0; m < 4; m++) {
                double b = bForMode(V, m);
                if (Double.isNaN(b)) System.out.printf("     -   ");
                else                  System.out.printf("   %.4f", b);
            }
            System.out.println();
        }
    }
}
