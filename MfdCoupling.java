/**
 * Mode field diameter (MFD) of a channel waveguide and butt-coupling loss to SMF.
 *
 * Reference mode is the SMF fundamental (circular Gaussian, w_smf = MFD_smf / 2).
 *
 *  - Cladding index from Sellmeier (Malitson 1965); core = clad / (1 - delta/100).
 *  - Waveguide mode is separable (EIM, TE-TE); each 1D profile is the symmetric-slab
 *    fundamental: cos(kappa*r) inside, cos(kappa*a)*exp(-gamma*(|r|-a)) outside.
 *  - MFD per direction (D4sigma):  sigma^2 = int r^2 I dr / int I dr,  MFD = 4*sigma.
 *  - Coupling loss to SMF:  eta = eta_x*eta_y,
 *       eta_dir = (int G*F)^2 / (int G^2 * int F^2),  loss_dB = -10 log10(eta).
 *
 * Compile: javac MfdCoupling.java
 * Run:     java MfdCoupling <wl_um> <MFD_smf_um> <W_um> <H_um> <delta_percent>
 *          java MfdCoupling 1.55 10.4 5.0 5.0 0.75
 */
public class MfdCoupling {

    public static double nSilica(double wl) {
        double B1 = 0.6961663, B2 = 0.4079426, B3 = 0.8974794;
        double C1 = 0.0684043*0.0684043, C2 = 0.1162414*0.1162414, C3 = 9.896161*9.896161;
        double l2 = wl*wl;
        return Math.sqrt(1 + B1*l2/(l2-C1) + B2*l2/(l2-C2) + B3*l2/(l2-C3));
    }

    /** {kappa, gamma, nEff} for symmetric-slab fundamental TE mode, or null. */
    public static double[] slabTeFundamental(double wl, double d, double nCore, double nClad) {
        if (nCore <= nClad) return null;
        double k0 = 2*Math.PI/wl, a = d/2.0;
        double V = k0 * a * Math.sqrt(nCore*nCore - nClad*nClad);
        if (V <= 0) return null;
        double upper = Math.min(Math.PI/2 - 1e-9, V - 1e-12);
        if (upper <= 0) return null;
        double lo = 1e-9, hi = upper, flo = fEig(lo, V), fhi = fEig(hi, V);
        if (flo*fhi > 0) return null;
        for (int i = 0; i < 100; i++) {
            double mid = 0.5*(lo+hi), fmid = fEig(mid, V);
            if (flo*fmid <= 0) { hi = mid; fhi = fmid; } else { lo = mid; flo = fmid; }
            if (hi - lo < 1e-13) break;
        }
        double u = 0.5*(lo+hi);
        double wv = Math.sqrt(Math.max(V*V - u*u, 0));
        double kappa = u/a, gamma = wv/a;
        double nEff = Math.sqrt(nCore*nCore - (kappa/k0)*(kappa/k0));
        return new double[] { kappa, gamma, nEff };
    }

    private static double fEig(double u, double V) {
        return u*Math.tan(u) - Math.sqrt(Math.max(V*V - u*u, 0));
    }

    private static double slabField(double r, double kappa, double gamma, double a) {
        return Math.abs(r) <= a ? Math.cos(kappa*r)
                                : Math.cos(kappa*a)*Math.exp(-gamma*(Math.abs(r)-a));
    }

    private static double gaussian(double r, double w) { return Math.exp(-(r*r)/(w*w)); }

    interface F1 { double at(double r); }

    private static double integrate(F1 f, double lo, double hi, int n) {
        double h = (hi-lo)/n, s = 0.5*(f.at(lo)+f.at(hi));
        for (int i = 1; i < n; i++) s += f.at(lo + i*h);
        return s*h;
    }

    private static double mfd1d(F1 field, double L, int N) {
        double num = integrate(r -> r*r*field.at(r)*field.at(r), -L, L, N);
        double den = integrate(r -> field.at(r)*field.at(r), -L, L, N);
        return den > 0 ? 4*Math.sqrt(num/den) : 0;
    }

    public static void main(String[] args) {
        if (args.length != 5) {
            System.out.println("Usage: java MfdCoupling <wl_um> <MFD_smf_um> <W_um> <H_um> <delta_percent>");
            System.out.println("Example: java MfdCoupling 1.55 10.4 5.0 5.0 0.75");
            return;
        }
        double wl = Double.parseDouble(args[0]);
        double mfd = Double.parseDouble(args[1]);
        double W = Double.parseDouble(args[2]);
        double H = Double.parseDouble(args[3]);
        double delta = Double.parseDouble(args[4]);

        double nClad = nSilica(wl);
        double nCore = nClad / (1 - delta/100);

        double[] vm = slabTeFundamental(wl, H, nCore, nClad);
        if (vm == null) { System.out.println("No guided mode under given parameters."); return; }
        double[] hm = slabTeFundamental(wl, W, vm[2], nClad);
        if (hm == null) { System.out.println("No guided mode under given parameters."); return; }

        final double kx = hm[0], gx = hm[1], ax = W/2;
        final double ky = vm[0], gy = vm[1], ay = H/2;
        final double wSmf = mfd/2;
        double L = Math.max(6*wSmf, Math.max(6*W, 6*H));
        int N = 4000;

        F1 fx = r -> slabField(r, kx, gx, ax);
        F1 fy = r -> slabField(r, ky, gy, ay);

        double mfdX = mfd1d(fx, L, N), mfdY = mfd1d(fy, L, N);
        double mfdEff = Math.sqrt(mfdX*mfdY);

        double etaX = overlap(fx, wSmf, L, N);
        double etaY = overlap(fy, wSmf, L, N);
        double eta = etaX*etaY;
        double lossDb = eta > 0 ? -10*Math.log10(eta) : Double.POSITIVE_INFINITY;

        System.out.printf("wavelength             = %.4f um%n", wl);
        System.out.printf("SMF MFD (reference)    = %.4f um%n", mfd);
        System.out.printf("WG core W x H          = %.4f x %.4f um%n", W, H);
        System.out.printf("index contrast         = %.4f %%%n", delta);
        System.out.printf("n_clad (Sellmeier)     = %.6f%n", nClad);
        System.out.printf("n_core (from delta)    = %.6f%n", nCore);
        System.out.printf("n_eff (vertical slab)  = %.6f%n", vm[2]);
        System.out.printf("WG MFD x (horizontal)  = %.4f um%n", mfdX);
        System.out.printf("WG MFD y (vertical)    = %.4f um%n", mfdY);
        System.out.printf("WG MFD (equiv circular)= %.4f um%n", mfdEff);
        System.out.printf("eta (total)            = %.6f%n", eta);
        System.out.printf("coupling loss (vs SMF) = %.3f dB%n", lossDb);
    }

    private static double overlap(F1 field, double wSmf, double L, int N) {
        double num = integrate(r -> gaussian(r, wSmf)*field.at(r), -L, L, N);
        double dg = integrate(r -> gaussian(r, wSmf)*gaussian(r, wSmf), -L, L, N);
        double df = integrate(r -> field.at(r)*field.at(r), -L, L, N);
        return (dg > 0 && df > 0) ? (num*num)/(dg*df) : 0;
    }
}
