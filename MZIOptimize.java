import java.util.*;
import java.util.function.ToDoubleFunction;

/**
 * MZI parameter optimizer (multi-start Nelder-Mead).
 * Finds (gap1, L1, gap2, L2, dL) minimizing sum of squared errors at target wavelengths.
 *
 *  Compile: javac MZIOptimize.java
 *  Run:     java MZIOptimize 5.0 5.0 0.75 0.98:1.0,1.55:0.0 [n_starts]
 */
public class MZIOptimize {

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
    public static double[] kappaNeff(double wl, double W, double H, double dl, double gap) {
        double nC = nSilica(wl), nCl = nC*(1-dl/100);
        double[] v = slabFundamentalTE(wl, H, nC, nCl); if (v == null) return null;
        double[] h = slabFundamentalTE(wl, W, v[2], nCl); if (h == null) return null;
        double k0 = 2*Math.PI/wl, beta = h[2]*k0, a = W/2, h2 = h[0]*h[0];
        double k = 2*h2*h[1]*Math.exp(-h[1]*gap) / (beta*(h2+h[1]*h[1])*(1+h[1]*a));
        return new double[] { k, h[2] };
    }
    public static Double pBar(double wl, double W, double H, double dl,
                              double g1, double L1, double g2, double L2, double dL) {
        double[] kn1 = kappaNeff(wl, W, H, dl, g1);
        double[] kn2 = kappaNeff(wl, W, H, dl, g2);
        if (kn1 == null || kn2 == null) return null;
        double kL1 = kn1[0]*L1, kL2 = kn2[0]*L2;
        double dphi = (2*Math.PI*kn1[1]/wl) * dL;
        double a_ = Math.cos(kL1)*Math.cos(kL2);
        double b_ = Math.sin(kL1)*Math.sin(kL2);
        return a_*a_ + b_*b_ - 2*a_*b_*Math.cos(dphi);
    }

    static double[][] BOUNDS = { {1, 15}, {10, 5000}, {1, 15}, {10, 5000}, {0.1, 200} };

    static double clip(double v, double lo, double hi) { return Math.max(lo, Math.min(hi, v)); }

    static double[] nelderMead(ToDoubleFunction<double[]> f, double[] x0, int maxIter) {
        int n = x0.length;
        final double[][] simp = new double[n+1][n];
        final double[] fv = new double[n+1];
        simp[0] = x0.clone();
        for (int i = 0; i < n; i++) {
            simp[i+1] = x0.clone();
            double step = Math.abs(x0[i]) > 1e-6 ? Math.abs(x0[i]) * 0.05 : 0.001;
            simp[i+1][i] += step;
        }
        for (int i = 0; i <= n; i++) fv[i] = f.applyAsDouble(simp[i]);

        Comparator<Integer> cmp = (a, b) -> Double.compare(fv[a], fv[b]);

        for (int it = 0; it < maxIter; it++) {
            Integer[] idx = new Integer[n+1];
            for (int i = 0; i <= n; i++) idx[i] = i;
            Arrays.sort(idx, cmp);
            double[][] ns = new double[n+1][];
            double[] nf = new double[n+1];
            for (int i = 0; i <= n; i++) { ns[i] = simp[idx[i]]; nf[i] = fv[idx[i]]; }
            for (int i = 0; i <= n; i++) { simp[i] = ns[i]; fv[i] = nf[i]; }
            if (fv[n] - fv[0] < 1e-8) break;
            double[] c = new double[n];
            for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) c[j] += simp[i][j] / n;
            double[] xr = new double[n];
            for (int j = 0; j < n; j++) xr[j] = c[j] + (c[j] - simp[n][j]);
            double fr = f.applyAsDouble(xr);
            if (fv[0] <= fr && fr < fv[n-1]) { simp[n] = xr; fv[n] = fr; continue; }
            if (fr < fv[0]) {
                double[] xe = new double[n];
                for (int j = 0; j < n; j++) xe[j] = c[j] + 2 * (xr[j] - c[j]);
                double fe = f.applyAsDouble(xe);
                if (fe < fr) { simp[n] = xe; fv[n] = fe; }
                else { simp[n] = xr; fv[n] = fr; }
                continue;
            }
            double[] xc = new double[n];
            for (int j = 0; j < n; j++) xc[j] = c[j] + 0.5 * (simp[n][j] - c[j]);
            double fc = f.applyAsDouble(xc);
            if (fc < fv[n]) { simp[n] = xc; fv[n] = fc; continue; }
            for (int i = 1; i <= n; i++) {
                for (int j = 0; j < n; j++) simp[i][j] = simp[0][j] + 0.5 * (simp[i][j] - simp[0][j]);
                fv[i] = f.applyAsDouble(simp[i]);
            }
        }
        Integer[] idx = new Integer[n+1];
        for (int i = 0; i <= n; i++) idx[i] = i;
        Arrays.sort(idx, cmp);
        return simp[idx[0]];
    }

    public static void main(String[] args) {
        if (args.length < 4) {
            System.out.println("Usage: java MZIOptimize <W> <H> <delta_pct> <targets> [n_starts=30]");
            System.out.println("Example: java MZIOptimize 5.0 5.0 0.75 0.98:1.0,1.55:0.0");
            return;
        }
        double W = Double.parseDouble(args[0]);
        double H = Double.parseDouble(args[1]);
        double dl = Double.parseDouble(args[2]);
        // Parse targets "wl:t[:w],..."
        List<double[]> targets = new ArrayList<>();
        for (String p : args[3].split(",")) {
            String[] b = p.split(":");
            double w = b.length >= 3 ? Double.parseDouble(b[2]) : 1.0;
            targets.add(new double[] { Double.parseDouble(b[0]), Double.parseDouble(b[1]), w });
        }
        int nStarts = args.length >= 5 ? Integer.parseInt(args[4]) : 30;

        Random rnd = new Random(1);
        ToDoubleFunction<double[]> cost = (x) -> {
            double[] xc = new double[x.length];
            double penalty = 0;
            for (int i = 0; i < x.length; i++) {
                xc[i] = clip(x[i], BOUNDS[i][0], BOUNDS[i][1]);
                double range = BOUNDS[i][1] - BOUNDS[i][0];
                penalty += 0.001 * Math.pow((x[i] - xc[i]) / Math.max(range, 1e-12), 2);
            }
            double s = 0;
            for (double[] t : targets) {
                Double p = pBar(t[0], W, H, dl, xc[0], xc[1], xc[2], xc[3], xc[4]);
                if (p == null) return 1e6;
                s += t[2] * (p - t[1]) * (p - t[1]);
            }
            return s + penalty;
        };

        double[] best = null;
        double bestC = Double.POSITIVE_INFINITY;
        for (int s = 0; s < nStarts; s++) {
            double[] x0 = new double[5];
            for (int i = 0; i < 5; i++) x0[i] = BOUNDS[i][0] + rnd.nextDouble() * (BOUNDS[i][1] - BOUNDS[i][0]);
            double[] x = nelderMead(cost, x0, 300);
            for (int i = 0; i < 5; i++) x[i] = clip(x[i], BOUNDS[i][0], BOUNDS[i][1]);
            double c = cost.applyAsDouble(x);
            if (c < bestC) { bestC = c; best = x; }
        }

        System.out.printf("best cost = %.6e%n", bestC);
        System.out.printf("DC1: gap = %.4f um   L = %.4f um%n", best[0], best[1]);
        System.out.printf("DC2: gap = %.4f um   L = %.4f um%n", best[2], best[3]);
        System.out.printf("Delta_L  = %.4f um%n%n", best[4]);
        System.out.println("verification (target -> achieved):");
        for (double[] t : targets) {
            Double p = pBar(t[0], W, H, dl, best[0], best[1], best[2], best[3], best[4]);
            System.out.printf("  lambda=%.4f um  P_bar target=%.3f  achieved=%.6f  P_cross=%.6f%n",
                              t[0], t[1], p, 1 - p);
        }
    }
}
