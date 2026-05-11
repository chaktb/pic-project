/**
 * WDM directional coupler — spectral table P_cross, P_through over wavelength.
 *
 *   kappa(lambda) varies with wavelength via Sellmeier n_core(lambda) and
 *   the slab modal parameters (gamma_x, beta).
 *   P_cross(lambda)   = sin^2(kappa(lambda) * L_DC)
 *   P_through(lambda) = cos^2(kappa(lambda) * L_DC)
 *
 * Compile: javac WDMCoupler.java
 * Run:     java WDMCoupler 1.50 1.60 5.0 5.0 0.75 5.0 600 50
 */
public class WDMCoupler {

    public static double nSilica(double wl) {
        double B1 = 0.6961663, B2 = 0.4079426, B3 = 0.8974794;
        double C1 = 0.0684043 * 0.0684043;
        double C2 = 0.1162414 * 0.1162414;
        double C3 = 9.896161 * 9.896161;
        double l2 = wl * wl;
        return Math.sqrt(1.0 + B1*l2/(l2-C1) + B2*l2/(l2-C2) + B3*l2/(l2-C3));
    }

    public static double[] slabFundamentalTE(double wl, double d, double nC, double nCl) {
        if (nC <= nCl) return null;
        double k0 = 2 * Math.PI / wl;
        double a = d / 2;
        double V = k0 * a * Math.sqrt(nC*nC - nCl*nCl);
        double upper = Math.min(Math.PI/2 - 1e-9, V - 1e-12);
        if (upper <= 0) return null;
        double lo = 1e-9, hi = upper, flo = f(lo, V), fhi = f(hi, V);
        if (flo*fhi > 0) return null;
        for (int i = 0; i < 100; i++) {
            double mid = 0.5*(lo+hi), fmid = f(mid, V);
            if (flo*fmid <= 0) { hi = mid; fhi = fmid; }
            else { lo = mid; flo = fmid; }
            if (hi - lo < 1e-13) break;
        }
        double u = 0.5*(lo+hi), w = Math.sqrt(Math.max(V*V - u*u, 0));
        return new double[] { u/a, w/a, Math.sqrt(nC*nC - (u/a/k0)*(u/a/k0)) };
    }
    private static double f(double u, double V) {
        return u * Math.tan(u) - Math.sqrt(Math.max(V*V - u*u, 0));
    }

    public static double kappaAt(double wl, double W, double H, double dl, double gap) {
        double nC = nSilica(wl);
        double nCl = nC * (1 - dl/100);
        double[] v = slabFundamentalTE(wl, H, nC, nCl);
        if (v == null) return Double.NaN;
        double[] h = slabFundamentalTE(wl, W, v[2], nCl);
        if (h == null) return Double.NaN;
        double kx = h[0], gx = h[1], nEff = h[2];
        double k0 = 2 * Math.PI / wl;
        double beta = nEff * k0;
        double a = W / 2;
        double h2 = kx * kx;
        return 2 * h2 * gx * Math.exp(-gx * gap) / (beta * (h2 + gx*gx) * (1 + gx * a));
    }

    public static void main(String[] args) {
        if (args.length < 7 || args.length > 8) {
            System.out.println("Usage: java WDMCoupler <wl_start> <wl_end> <W> <H> <delta_pct> <gap> <L_DC> [N=200]");
            return;
        }
        double wlA = Double.parseDouble(args[0]);
        double wlB = Double.parseDouble(args[1]);
        double W   = Double.parseDouble(args[2]);
        double H   = Double.parseDouble(args[3]);
        double dl  = Double.parseDouble(args[4]);
        double gap = Double.parseDouble(args[5]);
        double L   = Double.parseDouble(args[6]);
        int N      = args.length >= 8 ? Integer.parseInt(args[7]) : 200;

        System.out.printf("# WDM DC spectrum  W=%.2f H=%.2f delta=%.4f%% gap=%.2f L_DC=%.2f%n", W, H, dl, gap, L);
        System.out.println("# wavelength_um, P_cross, P_through");
        for (int i = 0; i < N; i++) {
            double wl = wlA + (wlB - wlA) * i / (N - 1);
            double k = kappaAt(wl, W, H, dl, gap);
            double pc = Double.isNaN(k) ? Double.NaN : Math.sin(k * L) * Math.sin(k * L);
            double pt = Double.isNaN(k) ? Double.NaN : Math.cos(k * L) * Math.cos(k * L);
            System.out.printf("%.6f, %.6f, %.6f%n", wl, pc, pt);
        }
    }
}
