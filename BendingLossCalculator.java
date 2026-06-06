/**
 * Bending (radiation) loss for a rectangular silica channel waveguide
 * (in-plane bend, Marcuse-style approximation, EIM-based mode parameters).
 *
 *  alpha(R) [1/um] = (gamma/2) * exp(-2*gamma*a) * exp(-B*R),
 *      B = 2*gamma^3 / (3*beta^2)
 *  Loss(dB) = 4.343 * alpha * R * theta
 *
 * Compile: javac BendingLossCalculator.java
 * Run:     java BendingLossCalculator 1.55 5.0 5.0 0.75 90 5
 */
public class BendingLossCalculator {

    public static double nSilica(double wl) {
        double B1 = 0.6961663, B2 = 0.4079426, B3 = 0.8974794;
        double C1 = 0.0684043 * 0.0684043;
        double C2 = 0.1162414 * 0.1162414;
        double C3 = 9.896161 * 9.896161;
        double l2 = wl * wl;
        return Math.sqrt(1.0 + B1*l2/(l2-C1) + B2*l2/(l2-C2) + B3*l2/(l2-C3));
    }

    /** Returns [kappa, gamma, n_eff] for fundamental TE; null if not guided. */
    public static double[] slabFundamentalTE(double wl, double d, double nCore, double nClad) {
        if (nCore <= nClad) return null;
        double k0 = 2 * Math.PI / wl;
        double a = d / 2.0;
        double V = k0 * a * Math.sqrt(nCore*nCore - nClad*nClad);
        if (V <= 0) return null;
        double upper = Math.min(Math.PI/2 - 1e-9, V - 1e-12);
        if (upper <= 0) return null;
        double lo = 1e-9, hi = upper;
        double flo = f(lo, V), fhi = f(hi, V);
        if (flo * fhi > 0) return null;
        for (int i = 0; i < 100; i++) {
            double mid = 0.5 * (lo + hi);
            double fmid = f(mid, V);
            if (flo * fmid <= 0) { hi = mid; fhi = fmid; }
            else { lo = mid; flo = fmid; }
            if (hi - lo < 1e-13) break;
        }
        double u = 0.5 * (lo + hi);
        double w = Math.sqrt(Math.max(V*V - u*u, 0.0));
        double kappa = u / a;
        double gamma = w / a;
        double nEff = Math.sqrt(nCore*nCore - (kappa/k0)*(kappa/k0));
        return new double[] { kappa, gamma, nEff };
    }

    private static double f(double u, double V) {
        return u * Math.tan(u) - Math.sqrt(Math.max(V*V - u*u, 0.0));
    }

    public static double bendingAlphaPerUm(double R, double gamma, double beta, double a) {
        double A = (gamma / 2.0) * Math.exp(-2.0 * gamma * a);
        double B = (2.0 * gamma * gamma * gamma) / (3.0 * beta * beta);
        return A * Math.exp(-B * R);
    }

    public static double bendingLossDb(double R, double theta, double gamma, double beta, double a) {
        return 4.342944819 * bendingAlphaPerUm(R, gamma, beta, a) * R * theta;
    }

    /** Smallest R [um] on the right-hand monotonic branch where loss <= target. */
    public static double criticalR(double targetDb, double theta, double gamma, double beta, double a) {
        double B = (2.0 * gamma * gamma * gamma) / (3.0 * beta * beta);
        if (B <= 0) return Double.NaN;
        double Rpeak = 1.0 / B;
        if (bendingLossDb(Rpeak, theta, gamma, beta, a) <= targetDb) return 0.0;
        double Rhi = Rpeak * 2;
        for (int i = 0; i < 80; i++) {
            if (bendingLossDb(Rhi, theta, gamma, beta, a) <= targetDb) break;
            Rhi *= 2;
            if (Rhi > 1e12) return Double.NaN;
        }
        double Rlo = Rpeak;
        for (int i = 0; i < 120; i++) {
            double mid = 0.5 * (Rlo + Rhi);
            if (bendingLossDb(mid, theta, gamma, beta, a) <= targetDb) Rhi = mid;
            else Rlo = mid;
            if (Rhi - Rlo < 1e-3) break;
        }
        return Rhi;
    }

    public static void main(String[] args) {
        if (args.length < 4 || args.length > 6) {
            System.out.println("Usage: java BendingLossCalculator <wl_um> <W_um> <H_um> <delta_pct> [theta_deg=90] [R_mm=5]");
            return;
        }
        double wl = Double.parseDouble(args[0]);
        double W  = Double.parseDouble(args[1]);
        double H  = Double.parseDouble(args[2]);
        double dl = Double.parseDouble(args[3]);
        double thetaDeg = args.length >= 5 ? Double.parseDouble(args[4]) : 90.0;
        double Rmm      = args.length >= 6 ? Double.parseDouble(args[5]) : 5.0;
        double theta = Math.toRadians(thetaDeg);
        double R = Rmm * 1000.0;

        double nClad = nSilica(wl);
        double nCore = nClad / (1.0 - dl/100.0);
        double[] vmode = slabFundamentalTE(wl, H, nCore, nClad);
        if (vmode == null) { System.out.println("No vertical mode."); return; }
        double[] hmode = slabFundamentalTE(wl, W, vmode[2], nClad);
        if (hmode == null) { System.out.println("No horizontal mode."); return; }
        double gamma = hmode[1], nEff = hmode[2];
        double k0 = 2 * Math.PI / wl;
        double beta = nEff * k0;
        double a = W / 2.0;

        double alpha = bendingAlphaPerUm(R, gamma, beta, a);
        double loss = bendingLossDb(R, theta, gamma, beta, a);
        double Rc01 = criticalR(0.1, theta, gamma, beta, a);
        double Rc1  = criticalR(1.0, theta, gamma, beta, a);

        System.out.printf("wavelength       = %.4f um%n", wl);
        System.out.printf("W x H            = %.4f x %.4f um%n", W, H);
        System.out.printf("index contrast   = %.4f %%%n", dl);
        System.out.printf("bend angle       = %.2f deg%n", thetaDeg);
        System.out.printf("R (point readout)= %.4f mm%n", Rmm);
        System.out.printf("n_core           = %.6f%n", nCore);
        System.out.printf("n_clad           = %.6f%n", nClad);
        System.out.printf("n_eff (total)    = %.6f%n", nEff);
        System.out.printf("gamma_x          = %.6f /um%n", gamma);
        System.out.printf("beta             = %.6f /um%n", beta);
        System.out.printf("alpha (at R)     = %.4e /um  (%.4e /mm)%n", alpha, alpha * 1e3);
        System.out.printf("loss @ R         = %.4f dB / %.0fdeg bend%n", loss, thetaDeg);
        System.out.printf("R_critical(0.1dB)= %.4f mm%n", Rc01 / 1000);
        System.out.printf("R_critical(1dB)  = %.4f mm%n", Rc1 / 1000);
    }
}
