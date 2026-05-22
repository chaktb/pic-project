/**
 * Bandpass filter from an asymmetric Mach-Zehnder interferometer.
 *
 * Two identical directional couplers (gap, L_DC) separated by arms with a
 * path-length difference Delta_L. The cross port forms a bandpass centered
 * where the arm phase Delta_phi = 0.
 *
 * Design from a filter spec (waveguide W/H/Delta + DC gap given):
 *   FSR      -> Delta_L = lambda_c^2 / (n_eff * FSR), snapped to an integer
 *               number of wavelengths so a cross peak sits at lambda_c.
 *   lambda_c -> peak location and the wavelength axis.
 *   passband -> coupler angle theta = kappa(lambda_c) * L_DC, from
 *                 P_cross = sin^2(2 theta) * cos^2(Delta_phi/2);
 *               the absolute -3 dB width spans 0..FSR/2 by under-coupling,
 *               trading peak transmission sin^2(2 theta) below FSR/2.
 *
 * Transfer (matches band_pass_filter.py / mzi):
 *   P_bar = a^2 + b^2 - 2ab cos(Delta_phi), a = cos(kL1)cos(kL2),
 *                                           b = sin(kL1)sin(kL2)
 *
 * Compile: javac BandPassFilter.java
 * Run:     java BandPassFilter <W> <H> <delta%> <gap> <lam_c_um> <FSR_nm> <PB_nm> [span_xFSR=3] [N=600]
 * Example: java BandPassFilter 5.0 5.0 0.75 5.0 1.55 1.6 0.6
 */
public class BandPassFilter {

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

    /** Bar-port power of the 2-coupler MZI, or NaN. */
    public static double mziPbar(double wl, double W, double H, double dl,
                                 double gap1, double L1, double gap2, double L2, double dL) {
        double[] kn1 = kappaNeff(wl, W, H, dl, gap1);
        double[] kn2 = kappaNeff(wl, W, H, dl, gap2);
        if (kn1 == null || kn2 == null) return Double.NaN;
        double kL1 = kn1[0] * L1, kL2 = kn2[0] * L2;
        double dphi = (2*Math.PI*kn1[1]/wl) * dL;
        double a_ = Math.cos(kL1) * Math.cos(kL2);
        double b_ = Math.sin(kL1) * Math.sin(kL2);
        return a_*a_ + b_*b_ - 2*a_*b_*Math.cos(dphi);
    }

    /** Returns [kappaC, neC, m, dL, fsrAchUm, theta, Ldc, peakT, r, clamped] or null. */
    public static double[] design(double W, double H, double dl, double gap,
                                  double lamC, double fsrNm, double pbNm) {
        double[] kn = kappaNeff(lamC, W, H, dl, gap);
        if (kn == null) return null;
        double kappaC = kn[0], neC = kn[1];
        double fsrUm = fsrNm / 1000.0, pbUm = pbNm / 1000.0;
        double dLtarget = lamC*lamC / (neC * fsrUm);
        double m = Math.max(1, Math.round(neC * dLtarget / lamC));
        double dL = m * lamC / neC;
        double fsrAchUm = lamC*lamC / (neC * dL);
        double rRaw = pbUm / fsrAchUm;
        double r = Math.min(Math.max(rRaw, 1e-6), 0.5);
        double theta = 0.5 * Math.asin(Math.min(1.0, Math.sqrt(0.5) / Math.cos(Math.PI * r / 2)));
        double Ldc = theta / kappaC;
        double peakT = Math.sin(2*theta) * Math.sin(2*theta);
        return new double[] { kappaC, neC, m, dL, fsrAchUm, theta, Ldc, peakT, r, rRaw > 0.5 ? 1 : 0 };
    }

    public static void main(String[] args) {
        if (args.length < 7) {
            System.out.println("Usage: java BandPassFilter <W> <H> <delta%> <gap>"
                + " <lam_c_um> <FSR_nm> <PB_nm> [span_xFSR=3] [N=600]");
            System.out.println("Example: java BandPassFilter 5.0 5.0 0.75 5.0 1.55 1.6 0.6");
            return;
        }
        double W = Double.parseDouble(args[0]);
        double H = Double.parseDouble(args[1]);
        double dl = Double.parseDouble(args[2]);
        double gap = Double.parseDouble(args[3]);
        double lamC = Double.parseDouble(args[4]);
        double fsrNm = Double.parseDouble(args[5]);
        double pbNm = Double.parseDouble(args[6]);
        double span = args.length >= 8 ? Double.parseDouble(args[7]) : 3.0;
        int N = args.length >= 9 ? Integer.parseInt(args[8]) : 600;

        double[] d = design(W, H, dl, gap, lamC, fsrNm, pbNm);
        if (d == null) {
            System.out.println("# design failed: no guided mode for this waveguide/wavelength");
            return;
        }
        double kappaC = d[0], neC = d[1];
        int m = (int) d[2];
        double dL = d[3], fsrAchUm = d[4], theta = d[5], Ldc = d[6], peakT = d[7];
        boolean clamped = d[9] != 0;

        System.out.printf("# Bandpass (asymmetric MZI)  W=%s H=%s delta=%s%% gap=%s%n", W, H, dl, gap);
        System.out.printf("# spec: lambda_c=%s um  FSR=%s nm  passband=%s nm%n", lamC, fsrNm, pbNm);
        System.out.printf("#   n_eff(lam_c)   = %.6f%n", neC);
        System.out.printf("#   kappa(lam_c)   = %.6e /um%n", kappaC);
        System.out.printf("#   Delta_L        = %.4f um   (m = %d wavelengths)%n", dL, m);
        System.out.printf("#   FSR achieved   = %.4f nm%n", fsrAchUm * 1000.0);
        System.out.printf("#   coupler angle  = %.4f rad (%.2f deg)%n", theta, Math.toDegrees(theta));
        System.out.printf("#   L_DC (x2)      = %.4f um   gap = %s um%n", Ldc, gap);
        System.out.printf("#   peak T (cross) = %.4f%n", peakT);
        if (clamped) {
            System.out.println("#   NOTE: requested passband > FSR/2; clamped to FSR/2 (3 dB couplers)."
                + " Use a cascade/lattice for narrower flat-top passbands.");
        }
        System.out.println("# wavelength_um, P_bar, P_cross");

        double wlA = lamC - 0.5 * span * fsrAchUm;
        double wlB = lamC + 0.5 * span * fsrAchUm;
        for (int i = 0; i < N; i++) {
            double wl = wlA + (wlB - wlA) * i / (N - 1);
            double pb = mziPbar(wl, W, H, dl, gap, Ldc, gap, Ldc, dL);
            if (Double.isNaN(pb)) System.out.printf("%.6f, NaN, NaN%n", wl);
            else System.out.printf("%.6f, %.6f, %.6f%n", wl, pb, Math.max(0.0, 1.0 - pb));
        }
    }
}
