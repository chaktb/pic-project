/**
 * Directional coupler (two identical parallel waveguides) — coupling and power transfer.
 *
 *   kappa = 2 * h^2 * q * exp(-q*gap) / ( beta * (h^2 + q^2) * (1 + q*a) )
 *   L_c   = pi / (2*kappa),   L_50 = pi / (4*kappa)
 *   P_cross(z)   = sin^2(kappa*z)
 *   P_through(z) = cos^2(kappa*z)
 *
 * Compile: javac DirectionalCoupler.java
 * Run:     java DirectionalCoupler 1.55 5.0 5.0 0.75 3.0 200
 */
public class DirectionalCoupler {

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

    public static void main(String[] args) {
        if (args.length != 6) {
            System.out.println("Usage: java DirectionalCoupler <wl_um> <W_um> <H_um> <delta_pct> <gap_um> <L_DC_um>");
            return;
        }
        double wl  = Double.parseDouble(args[0]);
        double W   = Double.parseDouble(args[1]);
        double H   = Double.parseDouble(args[2]);
        double dl  = Double.parseDouble(args[3]);
        double gap = Double.parseDouble(args[4]);
        double L_DC = Double.parseDouble(args[5]);

        double nCore = nSilica(wl);
        double nClad = nCore * (1.0 - dl / 100.0);
        double[] vmode = slabFundamentalTE(wl, H, nCore, nClad);
        if (vmode == null) { System.out.println("No vertical mode."); return; }
        double[] hmode = slabFundamentalTE(wl, W, vmode[2], nClad);
        if (hmode == null) { System.out.println("No horizontal mode."); return; }
        double kx = hmode[0], gx = hmode[1], nEff = hmode[2];
        double k0 = 2 * Math.PI / wl;
        double beta = nEff * k0;
        double a = W / 2;
        double h2 = kx * kx;
        double kappa = 2 * h2 * gx * Math.exp(-gx * gap) / (beta * (h2 + gx*gx) * (1 + gx * a));

        double Lc  = Math.PI / (2 * kappa);
        double L50 = Math.PI / (4 * kappa);
        double Pcross   = Math.pow(Math.sin(kappa * L_DC), 2);
        double Pthrough = Math.pow(Math.cos(kappa * L_DC), 2);

        System.out.printf("wavelength       = %.4f um%n", wl);
        System.out.printf("WG W x H         = %.4f x %.4f um%n", W, H);
        System.out.printf("index contrast   = %.4f %%%n", dl);
        System.out.printf("gap              = %.4f um%n", gap);
        System.out.printf("L_DC             = %.4f um   (%.4f mm)%n", L_DC, L_DC/1000);
        System.out.printf("n_core           = %.6f%n", nCore);
        System.out.printf("n_clad           = %.6f%n", nClad);
        System.out.printf("n_eff (total)    = %.6f%n", nEff);
        System.out.printf("kappa_x          = %.6f /um%n", kx);
        System.out.printf("gamma_x          = %.6f /um%n", gx);
        System.out.printf("beta             = %.6f /um%n", beta);
        System.out.printf("kappa (coupling) = %.6e /um  (%.6f /mm)%n", kappa, kappa*1000);
        System.out.printf("L_c (full)       = %.4f um   (%.4f mm)%n", Lc, Lc/1000);
        System.out.printf("L_50 (3 dB)      = %.4f um   (%.4f mm)%n", L50, L50/1000);
        System.out.printf("P_cross   @ L_DC = %.6f   (%.3f dB)%n", Pcross,   -10*Math.log10(Math.max(Pcross, 1e-30)));
        System.out.printf("P_through @ L_DC = %.6f   (%.3f dB)%n", Pthrough, -10*Math.log10(Math.max(Pthrough, 1e-30)));
    }
}
