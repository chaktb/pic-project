import java.util.*;
import java.util.function.ToDoubleFunction;

/**
 * Lattice-form bandpass filter optimizer.
 * N-stage lattice = (N+1) directional couplers separated by N equal-length
 * delay lines (path difference dL on bottom arm). Free parameters:
 *   L_k (k = 0 .. N), dL.
 *
 *  Compile: javac LatticeFilter.java
 *  Run:     java LatticeFilter <W> <H> <delta%> <wl_c_um> <passband_um> <N_stages>
 *  Example: java LatticeFilter 5.0 5.0 0.75 1.55 0.0002 3
 */
public class LatticeFilter {

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
        double nCl = nSilica(wl), nC = nCl / (1-dl/100);
        double[] v = slabFundamentalTE(wl, H, nC, nCl); if (v == null) return null;
        double[] h = slabFundamentalTE(wl, W, v[2], nCl); if (h == null) return null;
        double k0 = 2*Math.PI/wl, beta = h[2]*k0, a = W/2, h2 = h[0]*h[0];
        double k = 2*h2*h[1]*Math.exp(-h[1]*gap) / (beta*(h2+h[1]*h[1])*(1+h[1]*a));
        return new double[] { k, h[2] };
    }

    /** Bar power for the N-stage lattice. Ls length must be N+1. */
    public static Double latticePbar(double wl, double W, double H, double dl,
                                     double gap, double[] Ls, double dL) {
        double[] kn = kappaNeff(wl, W, H, dl, gap);
        if (kn == null) return null;
        double kappa = kn[0], nEff = kn[1];
        double dphi = (2*Math.PI*nEff/wl) * dL;
        double cp = Math.cos(dphi), sp = Math.sin(dphi);
        double uR=1, uI=0, vR=0, vI=0;
        int N = Ls.length - 1;
        for (int k = 0; k <= N; k++) {
            double theta = kappa * Ls[k];
            double c = Math.cos(theta), s = Math.sin(theta);
            double nuR = c*uR + s*vI;
            double nuI = c*uI - s*vR;
            double nvR = c*vR + s*uI;
            double nvI = c*vI - s*uR;
            uR = nuR; uI = nuI; vR = nvR; vI = nvI;
            if (k < N) {
                double nvR2 = vR*cp + vI*sp;
                double nvI2 = -vR*sp + vI*cp;
                vR = nvR2; vI = nvI2;
            }
        }
        return uR*uR + uI*uI;
    }

    static class Target { double wl, t, w; Target(double wl, double t, double w) {this.wl=wl; this.t=t; this.w=w;} }

    public static List<Target> bandpassTargets(double wlc, double pb) {
        return Arrays.asList(
            new Target(wlc,            1.00, 5.0),
            new Target(wlc - pb/2,     0.50, 1.0),
            new Target(wlc + pb/2,     0.50, 1.0),
            new Target(wlc - pb,       0.05, 2.0),
            new Target(wlc + pb,       0.05, 2.0),
            new Target(wlc - 2*pb,     0.00, 1.0),
            new Target(wlc + 2*pb,     0.00, 1.0)
        );
    }

    public static double[] clip(double[] x, double[][] b) {
        double[] r = new double[x.length];
        for (int i = 0; i < x.length; i++) r[i] = Math.max(b[i][0], Math.min(b[i][1], x[i]));
        return r;
    }

    public static double cost(double[] x, double W, double H, double dl, double gap,
                              int N, List<Target> targets, double[][] bounds) {
        double[] xc = clip(x, bounds);
        double[] Ls = Arrays.copyOf(xc, N+1);
        double dL = xc[N+1];
        double s = 0;
        for (Target t : targets) {
            Double p = latticePbar(t.wl, W, H, dl, gap, Ls, dL);
            if (p == null) return 1e6;
            s += t.w * (p - t.t) * (p - t.t);
        }
        for (int i = 0; i < x.length; i++) {
            double range = Math.max(bounds[i][1] - bounds[i][0], 1e-12);
            s += 0.001 * Math.pow((x[i] - xc[i]) / range, 2);
        }
        return s;
    }

    public static double[] nelderMead(ToDoubleFunction<double[]> f, double[] x0, int maxIter) {
        int n = x0.length;
        double[][] simp = new double[n+1][];
        simp[0] = x0.clone();
        for (int i = 0; i < n; i++) {
            simp[i+1] = x0.clone();
            double step = Math.abs(x0[i]) > 1e-6 ? Math.abs(x0[i])*0.05 : 0.001;
            simp[i+1][i] += step;
        }
        double[] fv = new double[n+1];
        for (int i = 0; i <= n; i++) fv[i] = f.applyAsDouble(simp[i]);
        for (int it = 0; it < maxIter; it++) {
            Integer[] idx = new Integer[n+1];
            for (int i = 0; i <= n; i++) idx[i] = i;
            final double[] fvSnap = fv;
            Arrays.sort(idx, (a,b) -> Double.compare(fvSnap[a], fvSnap[b]));
            double[][] s2 = new double[n+1][];
            double[] fv2 = new double[n+1];
            for (int i = 0; i <= n; i++) { s2[i] = simp[idx[i]]; fv2[i] = fv[idx[i]]; }
            simp = s2; fv = fv2;
            if (fv[n] - fv[0] < 1e-10) break;
            double[] c = new double[n];
            for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) c[j] += simp[i][j]/n;
            double[] xr = new double[n];
            for (int j = 0; j < n; j++) xr[j] = c[j] + (c[j] - simp[n][j]);
            double fr = f.applyAsDouble(xr);
            if (fv[0] <= fr && fr < fv[n-1]) { simp[n] = xr; fv[n] = fr; continue; }
            if (fr < fv[0]) {
                double[] xe = new double[n];
                for (int j = 0; j < n; j++) xe[j] = c[j] + 2*(xr[j] - c[j]);
                double fe = f.applyAsDouble(xe);
                if (fe < fr) { simp[n] = xe; fv[n] = fe; } else { simp[n] = xr; fv[n] = fr; }
                continue;
            }
            double[] xc = new double[n];
            for (int j = 0; j < n; j++) xc[j] = c[j] + 0.5*(simp[n][j] - c[j]);
            double fc = f.applyAsDouble(xc);
            if (fc < fv[n]) { simp[n] = xc; fv[n] = fc; continue; }
            for (int i = 1; i <= n; i++) {
                for (int j = 0; j < n; j++) simp[i][j] = simp[0][j] + 0.5*(simp[i][j] - simp[0][j]);
                fv[i] = f.applyAsDouble(simp[i]);
            }
        }
        int best = 0;
        for (int i = 1; i <= n; i++) if (fv[i] < fv[best]) best = i;
        return simp[best];
    }

    public static double[] optimize(double W, double H, double dl, double gap,
                                    int N, List<Target> targets,
                                    double[][] bounds, int nStarts, long seed) {
        Random rng = new Random(seed);
        ToDoubleFunction<double[]> f = x -> cost(x, W, H, dl, gap, N, targets, bounds);
        double[] best = null;
        double bestC = Double.POSITIVE_INFINITY;
        for (int s = 0; s < nStarts; s++) {
            double[] x0 = new double[bounds.length];
            for (int i = 0; i < bounds.length; i++) {
                x0[i] = bounds[i][0] + rng.nextDouble()*(bounds[i][1] - bounds[i][0]);
            }
            double[] r = nelderMead(f, x0, 400);
            r = clip(r, bounds);
            double c = f.applyAsDouble(r);
            if (c < bestC) { bestC = c; best = r; }
        }
        return best;
    }

    public static void main(String[] args) {
        if (args.length < 6) {
            System.out.println("Usage: java LatticeFilter <W> <H> <delta%> <wl_c_um> <passband_um> <N_stages>");
            return;
        }
        double W   = Double.parseDouble(args[0]);
        double H   = Double.parseDouble(args[1]);
        double dl  = Double.parseDouble(args[2]);
        double wlc = Double.parseDouble(args[3]);
        double pb  = Double.parseDouble(args[4]);
        int    N   = Integer.parseInt(args[5]);
        double gap = 5.0;

        List<Target> targets = bandpassTargets(wlc, pb);
        double[][] bounds = new double[N+2][2];
        for (int i = 0; i <= N; i++) { bounds[i][0] = 10; bounds[i][1] = 5000; }
        bounds[N+1][0] = 50; bounds[N+1][1] = 5000;
        double[] x = optimize(W, H, dl, gap, N, targets, bounds, 80, 1L);
        double[] Ls = Arrays.copyOf(x, N+1);
        double dL = x[N+1];

        System.out.printf("N stages = %d, center = %.4f um, passband = %.4f nm%n", N, wlc, pb*1000);
        for (int i = 0; i <= N; i++) System.out.printf("  DC%d: gap=%.3f um  L=%.4f um%n", i, gap, Ls[i]);
        System.out.printf("  dL  = %.4f um%n", dL);
        for (double off : new double[]{-2, -1, -0.5, 0, 0.5, 1, 2}) {
            double wl = wlc + off * pb;
            Double p = latticePbar(wl, W, H, dl, gap, Ls, dL);
            System.out.printf("  lambda=%.4f nm (%+.1f Δλ)  P_bar=%.6f%n", wl*1000, off, p);
        }
    }
}
