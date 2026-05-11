/**
 * MZI (Mach-Zehnder Interferometer) spectrum.
 * Two identical DCs (gap, L_DC) with arm path difference Delta_L.
 *
 *   E_bar   = cos(k1)cos(k2) - sin(k1)sin(k2) * exp(-j dphi)
 *   P_bar   = cos^2 + sin^2 - 2 cos sin cos(dphi) products as derived;
 *             code uses |E_bar|^2 = a^2 + b^2 - 2ab cos(dphi), a=cos(k1)cos(k2),
 *             b=sin(k1)sin(k2). P_cross = 1 - P_bar.
 *
 * Compile: javac MZI.java
 * Run:     java MZI 1.50 1.60 5.0 5.0 0.75 5.0 300 100 200
 */
public class MZI {

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

    /** Returns [kappa, n_eff] at wl, or null. */
    public static double[] kappaNeff(double wl, double W, double H, double dl, double gap) {
        double nC = nSilica(wl), nCl = nC * (1 - dl/100);
        double[] v = slabFundamentalTE(wl, H, nC, nCl);
        if (v == null) return null;
        double[] h = slabFundamentalTE(wl, W, v[2], nCl);
        if (h == null) return null;
        double kx = h[0], gx = h[1], nEff = h[2];
        double k0 = 2*Math.PI/wl, beta = nEff*k0, a = W/2, h2 = kx*kx;
        double kappa = 2*h2*gx*Math.exp(-gx*gap) / (beta*(h2+gx*gx)*(1+gx*a));
        return new double[] { kappa, nEff };
    }

    public static double[] mziPower(double wl, double W, double H, double dl, double gap, double L_DC, double dL) {
        double[] kn = kappaNeff(wl, W, H, dl, gap);
        if (kn == null) return null;
        double kappa = kn[0], nEff = kn[1];
        double kL = kappa * L_DC;
        double dphi = (2 * Math.PI * nEff / wl) * dL;
        double a_ = Math.cos(kL) * Math.cos(kL);
        double b_ = Math.sin(kL) * Math.sin(kL);
        double Pbar = a_*a_ + b_*b_ - 2*a_*b_*Math.cos(dphi);
        // Note: a_=cos^2(kL), b_=sin^2(kL); the cross terms in |E|^2 give:
        // |E_bar|^2 = (cos*cos)^2 + (sin*sin)^2 - 2*(cos*cos)*(sin*sin)*cos(dphi)
        // which is what's coded above with a_=cos*cos, b_=sin*sin
        return new double[] { kappa, nEff, kL, dphi, Pbar, Math.max(0, 1 - Pbar) };
    }

    public static void main(String[] args) {
        if (args.length < 8 || args.length > 9) {
            System.out.println("Usage: java MZI <wl_start> <wl_end> <W> <H> <delta_pct> <gap> <L_DC> <Delta_L> [N=400]");
            return;
        }
        double wlA = Double.parseDouble(args[0]);
        double wlB = Double.parseDouble(args[1]);
        double W   = Double.parseDouble(args[2]);
        double H   = Double.parseDouble(args[3]);
        double dl  = Double.parseDouble(args[4]);
        double gap = Double.parseDouble(args[5]);
        double L_DC = Double.parseDouble(args[6]);
        double dL  = Double.parseDouble(args[7]);
        int N      = args.length >= 9 ? Integer.parseInt(args[8]) : 400;

        double wlc = 0.5 * (wlA + wlB);
        double[] rc = mziPower(wlc, W, H, dl, gap, L_DC, dL);
        if (rc != null) {
            double fsr_um = wlc * wlc / (rc[1] * Math.max(dL, 1e-12));
            System.out.printf("# MZI W=%.2f H=%.2f delta=%.4f%% gap=%.2f L_DC=%.2f dL=%.2f%n", W, H, dl, gap, L_DC, dL);
            System.out.printf("# at lambda_center=%.4f: kappa=%.4e/um  kL=%.4f  n_eff=%.6f%n", wlc, rc[0], rc[2], rc[1]);
            System.out.printf("# approximate FSR ~ %.4f nm%n", fsr_um * 1000);
        }
        System.out.println("# wavelength_um, P_bar, P_cross");
        for (int i = 0; i < N; i++) {
            double wl = wlA + (wlB - wlA) * i / (N - 1);
            double[] r = mziPower(wl, W, H, dl, gap, L_DC, dL);
            if (r == null) System.out.printf("%.6f, NaN, NaN%n", wl);
            else           System.out.printf("%.6f, %.6f, %.6f%n", wl, r[4], r[5]);
        }
    }
}
