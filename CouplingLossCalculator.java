/**
 * Butt-coupling loss between SMF and a rectangular channel waveguide.
 * SMF: circular Gaussian (w = MFD/2). Waveguide: separable EIM mode.
 * Overlap integrals computed by trapezoidal rule on a fine grid.
 *
 * Compile: javac CouplingLossCalculator.java
 * Run:     java CouplingLossCalculator 1.55 10.4 5.0 5.0 0.75
 */
public class CouplingLossCalculator {

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

    public static double slabField(double r, double kappa, double gamma, double a) {
        if (Math.abs(r) <= a) return Math.cos(kappa * r);
        return Math.cos(kappa * a) * Math.exp(-gamma * (Math.abs(r) - a));
    }

    public static double gaussian(double r, double w) {
        return Math.exp(-(r * r) / (w * w));
    }

    /** Trapezoidal integration of f over [lo, hi] with N intervals. */
    public static double integrateTrap(java.util.function.DoubleUnaryOperator f,
                                       double lo, double hi, int N) {
        double h = (hi - lo) / N;
        double s = 0.5 * (f.applyAsDouble(lo) + f.applyAsDouble(hi));
        for (int i = 1; i < N; i++) s += f.applyAsDouble(lo + i * h);
        return s * h;
    }

    public static double[] couplingLoss(double wl, double mfd, double W, double H,
                                        double nCore, double nClad) {
        double[] vmode = slabFundamentalTE(wl, H, nCore, nClad);
        if (vmode == null) return null;
        double kappaY = vmode[0], gammaY = vmode[1], nEffV = vmode[2];
        double[] hmode = slabFundamentalTE(wl, W, nEffV, nClad);
        if (hmode == null) return null;
        double kappaX = hmode[0], gammaX = hmode[1];
        double aX = W / 2, aY = H / 2;
        double wSmf = mfd / 2;

        double L = Math.max(6 * wSmf, Math.max(6 * W, 6 * H));
        int N = 4000;

        // X direction
        double numX = integrateTrap(x -> gaussian(x, wSmf) * slabField(x, kappaX, gammaX, aX), -L, L, N);
        double dgX  = integrateTrap(x -> gaussian(x, wSmf) * gaussian(x, wSmf), -L, L, N);
        double dfX  = integrateTrap(x -> {
            double v = slabField(x, kappaX, gammaX, aX); return v * v;
        }, -L, L, N);
        double etaX = (numX * numX) / (dgX * dfX);

        // Y direction
        double numY = integrateTrap(y -> gaussian(y, wSmf) * slabField(y, kappaY, gammaY, aY), -L, L, N);
        double dgY  = integrateTrap(y -> gaussian(y, wSmf) * gaussian(y, wSmf), -L, L, N);
        double dfY  = integrateTrap(y -> {
            double v = slabField(y, kappaY, gammaY, aY); return v * v;
        }, -L, L, N);
        double etaY = (numY * numY) / (dgY * dfY);

        double eta = etaX * etaY;
        double lossDb = -10 * Math.log10(eta);
        return new double[] { nEffV, etaX, etaY, eta, lossDb };
    }

    public static void main(String[] args) {
        if (args.length != 5) {
            System.out.println("Usage: java CouplingLossCalculator <wavelength_um> <MFD_um> <W_um> <H_um> <delta_percent>");
            return;
        }
        double wl  = Double.parseDouble(args[0]);
        double mfd = Double.parseDouble(args[1]);
        double W   = Double.parseDouble(args[2]);
        double H   = Double.parseDouble(args[3]);
        double dl  = Double.parseDouble(args[4]);
        double nCore = nSilica(wl);
        double nClad = nCore * (1.0 - dl / 100.0);
        double[] r = couplingLoss(wl, mfd, W, H, nCore, nClad);
        if (r == null) { System.out.println("No guided mode."); return; }
        System.out.printf("wavelength            = %.4f um%n", wl);
        System.out.printf("SMF MFD               = %.4f um  (w_smf = %.4f um)%n", mfd, mfd/2);
        System.out.printf("WG core W x H         = %.4f x %.4f um%n", W, H);
        System.out.printf("index contrast        = %.4f %%%n", dl);
        System.out.printf("n_core (Sellmeier)    = %.6f%n", nCore);
        System.out.printf("n_clad                = %.6f%n", nClad);
        System.out.printf("n_eff (vertical slab) = %.6f%n", r[0]);
        System.out.printf("eta_x                 = %.6f%n", r[1]);
        System.out.printf("eta_y                 = %.6f%n", r[2]);
        System.out.printf("eta (total)           = %.6f%n", r[3]);
        System.out.printf("coupling loss         = %.3f dB%n", r[4]);
    }
}
