/**
 * Arrayed Waveguide Grating (AWG) - N x N wavelength demultiplexer.
 *
 * Two slab star couplers connected by an array of N_a waveguides with a
 * constant geometric path increment Delta_L. The grating order m fixes
 * Delta_L = m*lam_c/n_eff_a; the focal length f_l (Rowland circle) sets the
 * dispersion on the output face so adjacent channels are separated by
 * Delta_lam_ch.
 *
 *   Delta_L = m * lam_c / n_eff_a
 *   FSR_lam = lam_c^2 / (n_g * Delta_L)
 *   f_l     = (n_eff_s * d_aw * d_io * n_eff_a) / (m * n_g * Delta_lam_ch)
 *   theta_i = (i - (N_ch+1)/2) * d_io / f_l
 *   T_i(lam) = |sum_k A_k * exp(j*(k-c)*phi_i(lam))|^2 / (sum_k A_k)^2
 *   phi_i   = beta*Delta_L + 2*pi*n_eff_s*d_aw*theta_i/lam
 *   beta    = 2*pi*n_eff_a(lam) / lam
 *
 * Refractive indices: Sellmeier silica (Malitson 1965) + EIM (TE-TE).
 *
 * Compile: javac AWG.java
 * Run:     java AWG <W> <H> <delta%> <lam_c_um> <dlam_ch_nm> <N_ch> <N_a> <d_aw_um> <d_io_um>
 *                  [m=auto] [profile=gauss|uniform] [w_g=0.35] [span=1.0] [N=900]
 * Example: java AWG 6.0 6.0 0.75 1.55 0.8 4 30 15.0 20.0
 */
public class AWG {

    public static double nSilica(double wl) {
        double B1 = 0.6961663, B2 = 0.4079426, B3 = 0.8974794;
        double C1 = 0.0684043 * 0.0684043;
        double C2 = 0.1162414 * 0.1162414;
        double C3 = 9.896161 * 9.896161;
        double l2 = wl * wl;
        return Math.sqrt(1.0 + B1*l2/(l2-C1) + B2*l2/(l2-C2) + B3*l2/(l2-C3));
    }

    /** Returns the fundamental TE slab mode (kappa, n_eff), or null if cut off. */
    public static double[] slabFundamentalTE(double wl, double d, double nC, double nCl) {
        if (nC <= nCl) return null;
        double k0 = 2 * Math.PI / wl;
        double a = d / 2;
        double V = k0 * a * Math.sqrt(nC*nC - nCl*nCl);
        double upper = Math.min(Math.PI/2 - 1e-9, V - 1e-12);
        if (upper <= 0) return null;
        double lo = 1e-9, hi = upper;
        double flo = f(lo, V), fhi = f(hi, V);
        if (flo * fhi > 0) return null;
        for (int i = 0; i < 100; i++) {
            double mid = 0.5 * (lo + hi), fmid = f(mid, V);
            if (flo * fmid <= 0) { hi = mid; fhi = fmid; }
            else { lo = mid; flo = fmid; }
            if (hi - lo < 1e-13) break;
        }
        double u = 0.5 * (lo + hi);
        return new double[] { u/a, Math.sqrt(nC*nC - (u/a/k0)*(u/a/k0)) };
    }
    private static double f(double u, double V) {
        return u * Math.tan(u) - Math.sqrt(Math.max(V*V - u*u, 0));
    }

    /** EIM (TE-TE): returns [n_eff_a, n_eff_s] or null. */
    public static double[] eim(double wl, double W, double H, double deltaPct) {
        double nCore = nSilica(wl);
        double nClad = nCore * (1 - deltaPct / 100.0);
        double[] v = slabFundamentalTE(wl, H, nCore, nClad);
        if (v == null) return null;
        double nEffS = v[1];
        double[] h = slabFundamentalTE(wl, W, nEffS, nClad);
        if (h == null) return null;
        return new double[] { h[1], nEffS };
    }

    public static double groupIndex(double wl, double W, double H, double deltaPct) {
        double dwl = wl * 1e-4;
        double[] a = eim(wl - dwl, W, H, deltaPct);
        double[] b = eim(wl + dwl, W, H, deltaPct);
        double[] c = eim(wl, W, H, deltaPct);
        if (a == null || b == null || c == null) return Double.NaN;
        double dn = (b[0] - a[0]) / (2 * dwl);
        return c[0] - wl * dn;
    }

    public static class Design {
        double nEffA, nEffS, nG, dL, fsr, fl, sweep;
        int m;
    }

    public static Design design(double W, double H, double deltaPct, double lamC,
                                 double dlamChUm, int nCh, int nA,
                                 double dAw, double dIo, int mPref) {
        double[] at = eim(lamC, W, H, deltaPct);
        if (at == null) return null;
        double nEffA = at[0], nEffS = at[1];
        double nG = groupIndex(lamC, W, H, deltaPct);
        if (Double.isNaN(nG)) nG = nEffA;

        int m;
        if (mPref <= 0) {
            double fsrTarget = (nCh + 2) * dlamChUm;
            m = Math.max(1, (int) Math.round(lamC * nEffA / (nG * fsrTarget)));
        } else {
            m = mPref;
        }
        double dL = m * lamC / nEffA;
        double fsr = lamC * lamC / (nG * dL);
        double fl = (nEffS * dAw * dIo * nEffA) / (m * nG * dlamChUm);
        double sweep = (nCh - 1) * dIo / fl;

        Design d = new Design();
        d.nEffA = nEffA; d.nEffS = nEffS; d.nG = nG;
        d.m = m; d.dL = dL; d.fsr = fsr; d.fl = fl; d.sweep = sweep;
        return d;
    }

    public static double[] couplingAmp(int nA, String profile, double wG) {
        double[] A = new double[nA + 1];
        double c = 0.5 * (nA + 1);
        for (int k = 1; k <= nA; k++) {
            if ("gauss".equals(profile)) {
                double x = (k - c) / (wG * nA);
                A[k] = Math.exp(-x * x);
            } else {
                A[k] = 1.0;
            }
        }
        return A;
    }

    public static double transmission(double wl, int port, double W, double H, double deltaPct,
                                       double dAw, double dIo, int nA, int nCh,
                                       Design d, double[] A) {
        double[] at = eim(wl, W, H, deltaPct);
        if (at == null) return Double.NaN;
        double nEffA = at[0], nEffS = at[1];
        double beta = 2 * Math.PI * nEffA / wl;
        double cch = 0.5 * (nCh + 1);
        double theta = (port - cch) * dIo / d.fl;
        double dphi = beta * d.dL + 2 * Math.PI * nEffS * dAw * theta / wl;
        double c = 0.5 * (nA + 1);
        double re = 0, im = 0, sumA = 0;
        for (int k = 1; k <= nA; k++) {
            double phi = (k - c) * dphi;
            re += A[k] * Math.cos(phi);
            im += A[k] * Math.sin(phi);
            sumA += A[k];
        }
        return (re * re + im * im) / (sumA * sumA);
    }

    public static void main(String[] args) {
        if (args.length < 9) {
            System.out.println("Usage: java AWG <W> <H> <delta%> <lam_c_um> <dlam_ch_nm> <N_ch> <N_a>");
            System.out.println("              <d_aw_um> <d_io_um> [m=auto] [profile=gauss|uniform]");
            System.out.println("              [w_g=0.35] [span=1.0] [N=900]");
            return;
        }
        double W = Double.parseDouble(args[0]);
        double H = Double.parseDouble(args[1]);
        double deltaPct = Double.parseDouble(args[2]);
        double lamC = Double.parseDouble(args[3]);
        double dlamChNm = Double.parseDouble(args[4]);
        int nCh = Integer.parseInt(args[5]);
        int nA = Integer.parseInt(args[6]);
        double dAw = Double.parseDouble(args[7]);
        double dIo = Double.parseDouble(args[8]);
        int mPref = args.length > 9 ? Integer.parseInt(args[9]) : 0;
        String profile = args.length > 10 ? args[10] : "gauss";
        double wG = args.length > 11 ? Double.parseDouble(args[11]) : 0.35;
        double span = args.length > 12 ? Double.parseDouble(args[12]) : 1.0;
        int N = args.length > 13 ? Integer.parseInt(args[13]) : 900;

        double dlamChUm = dlamChNm * 1e-3;
        Design d = design(W, H, deltaPct, lamC, dlamChUm, nCh, nA, dAw, dIo, mPref);
        if (d == null) {
            System.out.println("No guided mode for these parameters.");
            return;
        }
        double[] A = couplingAmp(nA, profile, wG);

        System.out.printf("n_eff_a = %.6f%n", d.nEffA);
        System.out.printf("n_eff_s = %.6f%n", d.nEffS);
        System.out.printf("n_g     = %.6f%n", d.nG);
        System.out.printf("m       = %d%n", d.m);
        System.out.printf("Delta_L = %.4f um%n", d.dL);
        System.out.printf("FSR     = %.4f nm%n", d.fsr * 1000);
        System.out.printf("f_l     = %.4f mm%n", d.fl * 1e-3);
        System.out.printf("sweep   = %.3f deg%n", Math.toDegrees(d.sweep));

        System.out.println();
        StringBuilder header = new StringBuilder("# wavelength_um");
        for (int i = 1; i <= nCh; i++) header.append(String.format("  port%d_dB", i));
        System.out.println(header);
        double half = 0.5 * span * d.fsr;
        for (int n = 0; n < N; n++) {
            double wl = lamC - half + 2 * half * n / (N - 1);
            StringBuilder row = new StringBuilder(String.format("%.6f", wl));
            for (int i = 1; i <= nCh; i++) {
                double T = transmission(wl, i, W, H, deltaPct, dAw, dIo, nA, nCh, d, A);
                double dB = T > 1e-12 ? 10 * Math.log10(T) : -120.0;
                row.append(String.format("  %8.3f", dB));
            }
            System.out.println(row);
        }
    }
}
