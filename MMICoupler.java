/**
 * M x N MMI coupler dimension calculator (general interference, Soldano-Pennings).
 *
 *   pitch  = W_wg + gap
 *   W_MMI  = max(M, N) * pitch
 *   W_e    = W_MMI + (lambda/pi) * (n_clad/n_core)^(2*sigma) / sqrt(n_core^2 - n_clad^2)
 *   L_pi   = 4 * n_core * W_e^2 / (3 * lambda)
 *   L_MMI  = 3 * L_pi / N
 *
 * Compile: javac MMICoupler.java
 * Run:     java MMICoupler 1.55 5.0 5.0 0.75 2 2 5.0
 */
public class MMICoupler {

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
        if (args.length != 7) {
            System.out.println("Usage: java MMICoupler <wl_um> <W_wg_um> <H_um> <delta_pct> <M> <N> <gap_um>");
            return;
        }
        double wl   = Double.parseDouble(args[0]);
        double W_wg = Double.parseDouble(args[1]);
        double H    = Double.parseDouble(args[2]);
        double dl   = Double.parseDouble(args[3]);
        int M       = Integer.parseInt(args[4]);
        int N       = Integer.parseInt(args[5]);
        double gap  = Double.parseDouble(args[6]);
        double sigma = 0.0; // TE

        double nClad = nSilica(wl);
        double nCore = nClad / (1.0 - dl/100.0);
        double[] vmode = slabFundamentalTE(wl, H, nCore, nClad);
        double nEff = vmode != null ? vmode[2] : nCore;

        int Nmax = Math.max(M, N);
        double pitch = W_wg + gap;
        double W_MMI = Nmax * pitch;
        double W_e = W_MMI + (wl/Math.PI) * Math.pow(nClad/nCore, 2*sigma) / Math.sqrt(nCore*nCore - nClad*nClad);
        double L_pi = 4 * nCore * W_e * W_e / (3 * wl);
        double L_MMI = 3 * L_pi / N;

        System.out.printf("wavelength      = %.4f um%n", wl);
        System.out.printf("WG W x H        = %.4f x %.4f um%n", W_wg, H);
        System.out.printf("index contrast  = %.4f %%%n", dl);
        System.out.printf("M x N           = %d x %d%n", M, N);
        System.out.printf("gap             = %.4f um%n", gap);
        System.out.printf("n_core          = %.6f%n", nCore);
        System.out.printf("n_clad          = %.6f%n", nClad);
        System.out.printf("n_eff (vertical)= %.6f%n", nEff);
        System.out.printf("pitch           = %.4f um%n", pitch);
        System.out.printf("W_MMI           = %.4f um%n", W_MMI);
        System.out.printf("W_e (eff width) = %.4f um%n", W_e);
        System.out.printf("L_pi (beat)     = %.4f um   (%.4f mm)%n", L_pi, L_pi/1000);
        System.out.printf("L_MMI           = %.4f um   (%.4f mm)%n", L_MMI, L_MMI/1000);

        double y0 = -(Nmax - 1) / 2.0 * pitch;
        StringBuilder sb = new StringBuilder("[");
        for (int i = 0; i < M; i++) {
            sb.append(String.format("%.3f", y0 + (Nmax - M)/2.0*pitch + i*pitch));
            if (i < M-1) sb.append(", ");
        }
        sb.append("]");
        System.out.println("input  port y   = " + sb);
        sb = new StringBuilder("[");
        for (int i = 0; i < N; i++) {
            sb.append(String.format("%.3f", y0 + (Nmax - N)/2.0*pitch + i*pitch));
            if (i < N-1) sb.append(", ");
        }
        sb.append("]");
        System.out.println("output port y   = " + sb);
    }
}
