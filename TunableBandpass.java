/**
 * Thermally tunable bandpass filter (asymmetric MZI + NiCr thin-film heater).
 *
 * Two identical, tunable directional couplers (each a balanced-MZI coupler
 * under a NiCr heater) separated by arms with a fixed path-length difference
 * Delta_L. Delta_L is snapped to an integer number of wavelengths at lambda_c,
 * so the cross-port peak stays locked at lambda_c (center does NOT move). The
 * NiCr heater changes only the coupler angle theta, hence only the passband
 * WIDTH.
 *
 *   R = rho*L_h/(w_h*t_h);  P = V^2/R;  dT = R_th*P;  dn = (dn/dT)*dT
 *   d_theta = pi*dn*L_h/lambda_c;  theta(T) = theta_base + d_theta
 *   P_cross = sin^2(2 theta) * cos^2(Delta_phi/2)
 *   PB = FSR*(2/pi)*acos(sqrt(0.5/sin^2(2 theta)))    (0 if peak <= 0.5)
 *
 * Compile: javac TunableBandpass.java
 * Run:     java TunableBandpass <W> <H> <delta%> <gap> <lam_c_um> <FSR_nm> <base_PB_nm> <V> [span_xFSR=3] [N=600]
 * Example: java TunableBandpass 5.0 5.0 0.75 5.0 1.55 1.6 0.3 5.0
 */
public class TunableBandpass {

    // NiCr / silica defaults.
    static final double RHO_NICR = 1.1e-6;  // Ohm*m
    static final double DNDT = 1.0e-5;      // 1/K (fused silica)
    static final double R_TH = 2000.0;      // K/W
    static final double W_H = 4.0;          // um (heater width)
    static final double T_H = 0.1;          // um (film thickness)

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

    /** Bar-port power; dtheta is the heater-induced extra coupling angle. NaN if no mode. */
    public static double mziPbar(double wl, double W, double H, double dl,
                                 double gap1, double L1, double gap2, double L2,
                                 double dL, double dtheta) {
        double[] kn1 = kappaNeff(wl, W, H, dl, gap1);
        double[] kn2 = kappaNeff(wl, W, H, dl, gap2);
        if (kn1 == null || kn2 == null) return Double.NaN;
        double kL1 = kn1[0] * L1 + dtheta, kL2 = kn2[0] * L2 + dtheta;
        double dphi = (2*Math.PI*kn1[1]/wl) * dL;
        double a_ = Math.cos(kL1) * Math.cos(kL2);
        double b_ = Math.sin(kL1) * Math.sin(kL2);
        return a_*a_ + b_*b_ - 2*a_*b_*Math.cos(dphi);
    }

    /** Absolute -3 dB passband width (same units as fsr) for coupler angle theta. */
    public static double pbFromTheta(double theta, double fsr) {
        double peak = Math.sin(2*theta) * Math.sin(2*theta);
        if (peak <= 0.5) return 0.0;
        return fsr * (2.0/Math.PI) * Math.acos(Math.sqrt(0.5/peak));
    }

    public static void main(String[] args) {
        if (args.length < 8) {
            System.out.println("Usage: java TunableBandpass <W> <H> <delta%> <gap>"
                + " <lam_c_um> <FSR_nm> <base_PB_nm> <V> [span_xFSR=3] [N=600]");
            System.out.println("Example: java TunableBandpass 5.0 5.0 0.75 5.0 1.55 1.6 0.3 5.0");
            return;
        }
        double W = Double.parseDouble(args[0]);
        double H = Double.parseDouble(args[1]);
        double dl = Double.parseDouble(args[2]);
        double gap = Double.parseDouble(args[3]);
        double lamC = Double.parseDouble(args[4]);
        double fsrNm = Double.parseDouble(args[5]);
        double basePbNm = Double.parseDouble(args[6]);
        double volt = Double.parseDouble(args[7]);
        double span = args.length >= 9 ? Double.parseDouble(args[8]) : 3.0;
        int N = args.length >= 10 ? Integer.parseInt(args[9]) : 600;

        double[] kn = kappaNeff(lamC, W, H, dl, gap);
        if (kn == null) {
            System.out.println("# design failed: no guided mode for this waveguide/wavelength");
            return;
        }
        double kappaC = kn[0], neC = kn[1];
        double fsrUm = fsrNm / 1000.0, basePbUm = basePbNm / 1000.0;

        // Fixed center: snap Delta_L to an integer number of wavelengths at lamC.
        double dLtarget = lamC*lamC / (neC * fsrUm);
        long m = Math.max(1, Math.round(neC * dLtarget / lamC));
        double dL = m * lamC / neC;
        double fsrAchUm = lamC*lamC / (neC * dL);

        // Base coupler angle from the unpowered passband (<= FSR/2).
        double r0 = Math.min(Math.max(basePbUm / fsrAchUm, 1e-6), 0.5);
        double thetaBase = 0.5 * Math.asin(Math.min(1.0, Math.sqrt(0.5) / Math.cos(Math.PI * r0 / 2)));
        double Ldc = thetaBase / kappaC;

        // NiCr heater (length = L_DC).
        double lM = Ldc*1e-6, wM = W_H*1e-6, tM = T_H*1e-6;
        double R = RHO_NICR * lM / (wM * tM);
        double P = R > 0 ? volt*volt / R : 0.0;
        double dT = R_TH * P;
        double dn = DNDT * dT;
        double dtheta = Math.PI * dn * Ldc / lamC;
        double thetaT = thetaBase + dtheta;
        double peakT = Math.sin(2*thetaT) * Math.sin(2*thetaT);

        System.out.printf("# Tunable bandpass (asymmetric MZI + NiCr heater)  W=%s H=%s delta=%s%% gap=%s%n", W, H, dl, gap);
        System.out.printf("# spec: lambda_c=%s um (FIXED)  FSR=%s nm  base_PB=%s nm  V=%s V%n", lamC, fsrNm, basePbNm, volt);
        System.out.printf("#   n_eff(lam_c)   = %.6f%n", neC);
        System.out.printf("#   Delta_L        = %.4f um   (m = %d wavelengths, center locked)%n", dL, m);
        System.out.printf("#   FSR achieved   = %.4f nm%n", fsrAchUm * 1000.0);
        System.out.printf("#   L_DC (x2)      = %.4f um   (heater length = L_DC)%n", Ldc);
        System.out.printf("#   theta_base     = %.3f deg   base PB = %.4f nm%n",
                Math.toDegrees(thetaBase), pbFromTheta(thetaBase, fsrAchUm) * 1000.0);
        System.out.println("#   -- NiCr heater --");
        System.out.printf("#   R              = %.2f Ohm%n", R);
        System.out.printf("#   P = V^2/R      = %.4f mW%n", P * 1e3);
        System.out.printf("#   dT             = %.3f K%n", dT);
        System.out.printf("#   dn (T-O)       = %.3e%n", dn);
        System.out.printf("#   d_theta        = %.3f deg%n", Math.toDegrees(dtheta));
        System.out.printf("#   theta(T)       = %.3f deg%n", Math.toDegrees(thetaT));
        System.out.printf("#   peak T (cross) = %.4f%n", peakT);
        System.out.printf("#   PB tuned       = %.4f nm   (center unchanged at %s um)%n",
                pbFromTheta(thetaT, fsrAchUm) * 1000.0, lamC);
        System.out.println("# wavelength_um, P_bar, P_cross");

        double wlA = lamC - 0.5 * span * fsrAchUm;
        double wlB = lamC + 0.5 * span * fsrAchUm;
        for (int i = 0; i < N; i++) {
            double wl = wlA + (wlB - wlA) * i / (N - 1);
            double pb = mziPbar(wl, W, H, dl, gap, Ldc, gap, Ldc, dL, dtheta);
            if (Double.isNaN(pb)) System.out.printf("%.6f, NaN, NaN%n", wl);
            else System.out.printf("%.6f, %.6f, %.6f%n", wl, pb, Math.max(0.0, 1.0 - pb));
        }
    }
}
