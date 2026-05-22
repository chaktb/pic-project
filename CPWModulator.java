/**
 * LiNbO3 coplanar-waveguide (CPW) traveling-wave modulator calculator.
 *
 * Port of the physics in public/cpw_modulator.html. Over an RF sweep it
 * computes the CPW RLGC line parameters and S11/S21 (electrical bandwidth),
 * RF attenuation and complex characteristic impedance, microwave/optical
 * velocity (effective-index) matching, the traveling-wave electro-optic (EO)
 * response, and the half-wave voltage Vpi(f). The optical index combines
 * LiNbO3 Sellmeier dispersion with a Ti-indiffusion slab (fundamental mode);
 * quasi-static CPW parameters use the elliptic-integral conformal map.
 *
 *  Compile: javac CPWModulator.java
 *  Run:     java CPWModulator [key=value ...]
 *  Example: java CPWModulator w=8 s=16 line_len=20 pol=TM
 *  Keys (units): w,s,t_au,tox[um] line_len[mm] sigma_au[S/m] epsr_sio2 epsr_ln
 *                tan_sio2 tan_ln lam[um] dn_ti d_ti[um] pol(TM|TE) r33,r13[pm/V]
 *                gamma pp(0|1) f_start,f_stop[GHz] npts z0_port[ohm]
 */
import java.util.HashMap;
import java.util.Map;

public class CPWModulator {

    static final double MU0 = 4 * Math.PI * 1e-7;
    static final double EPS0 = 8.854187817e-12;
    static final double C0 = 1 / Math.sqrt(MU0 * EPS0);

    // -------------------- minimal complex --------------------
    static final class Cx {
        final double re, im;
        Cx(double re, double im) { this.re = re; this.im = im; }
        Cx add(Cx b) { return new Cx(re + b.re, im + b.im); }
        Cx sub(Cx b) { return new Cx(re - b.re, im - b.im); }
        Cx mul(Cx b) { return new Cx(re * b.re - im * b.im, re * b.im + im * b.re); }
        Cx div(Cx b) {
            double d = b.re * b.re + b.im * b.im;
            return new Cx((re * b.re + im * b.im) / d, (im * b.re - re * b.im) / d);
        }
        double abs() { return Math.hypot(re, im); }
        double arg() { return Math.atan2(im, re); }
        Cx sqrt() {
            double r = Math.sqrt(abs()), t = arg() / 2;
            return new Cx(r * Math.cos(t), r * Math.sin(t));
        }
        Cx exp() {
            double e = Math.exp(re);
            return new Cx(e * Math.cos(im), e * Math.sin(im));
        }
        Cx cosh() { return new Cx(Math.cosh(re) * Math.cos(im), Math.sinh(re) * Math.sin(im)); }
        Cx sinh() { return new Cx(Math.sinh(re) * Math.cos(im), Math.cosh(re) * Math.sin(im)); }
    }

    // -------------------- math utilities --------------------
    static double[] linspace(double a, double b, int n) {
        if (n <= 1) return new double[] { a };
        double s = (b - a) / (n - 1);
        double[] out = new double[n];
        for (int i = 0; i < n; i++) out[i] = a + i * s;
        return out;
    }
    static double interp(double x, double[] xs, double[] ys) {
        if (xs.length < 2) return ys.length > 0 ? ys[0] : 0;
        if (x <= xs[0]) return ys[0];
        if (x >= xs[xs.length - 1]) return ys[ys.length - 1];
        int lo = 0, hi = xs.length - 1;
        while (hi - lo > 1) {
            int m = (lo + hi) >> 1;
            if (xs[m] <= x) lo = m; else hi = m;
        }
        double t = (x - xs[lo]) / (xs[hi] - xs[lo]);
        return ys[lo] + t * (ys[hi] - ys[lo]);
    }
    /** First x where y crosses target; NaN if none. */
    static double findCrossing(double[] x, double[] y, double target) {
        for (int i = 0; i < x.length - 1; i++) {
            double d1 = y[i] - target, d2 = y[i + 1] - target;
            if (d1 == 0) return x[i];
            if (d1 * d2 < 0) {
                if (y[i + 1] == y[i]) return x[i];
                return x[i] + (target - y[i]) * (x[i + 1] - x[i]) / (y[i + 1] - y[i]);
            }
        }
        return Double.NaN;
    }
    static double ellipk(double m) {
        if (m >= 1) return Double.POSITIVE_INFINITY;
        if (m < 0) return Double.NaN;
        double a = 1, b = Math.sqrt(1 - m);
        for (int i = 0; i < 200; i++) {
            double an = 0.5 * (a + b), bn = Math.sqrt(a * b);
            if (Math.abs(an - bn) <= 1e-15 * an) { a = an; break; }
            a = an; b = bn;
        }
        return Math.PI / (2 * a);
    }

    interface Fn { double f(double x); }

    static Double brentq(Fn fn, double xa, double xb, double xtol, int maxIt) {
        double a = xa, b = xb, fa = fn.f(a), fb = fn.f(b);
        if (fa * fb > 0) return null;
        if (Math.abs(fa) < Math.abs(fb)) { double tmp = a; a = b; b = tmp; tmp = fa; fa = fb; fb = tmp; }
        double c = a, fc = fa, d = a;
        boolean mflag = true;
        for (int i = 0; i < maxIt; i++) {
            if (Math.abs(b - a) < xtol) return b;
            if (fb == 0) return b;
            double s;
            if (fa != fc && fb != fc) {
                s = a * fb * fc / ((fa - fb) * (fa - fc))
                  + b * fa * fc / ((fb - fa) * (fb - fc))
                  + c * fa * fb / ((fc - fa) * (fc - fb));
            } else {
                s = b - fb * (b - a) / (fb - fa);
            }
            double bnd = 0.25 * (3 * a + b);
            boolean c1 = (s - bnd) * (s - b) > 0;
            boolean c2 = mflag && Math.abs(s - b) >= 0.5 * Math.abs(b - c);
            boolean c3 = !mflag && Math.abs(s - b) >= 0.5 * Math.abs(c - d);
            boolean c4 = mflag && Math.abs(b - c) < xtol;
            boolean c5 = !mflag && Math.abs(c - d) < xtol;
            if (c1 || c2 || c3 || c4 || c5) { s = 0.5 * (a + b); mflag = true; }
            else mflag = false;
            double fs = fn.f(s);
            d = c; c = b; fc = fb;
            if (fa * fs < 0) { b = s; fb = fs; } else { a = s; fa = fs; }
            if (Math.abs(fa) < Math.abs(fb)) { double tmp = a; a = b; b = tmp; tmp = fa; fa = fb; fb = tmp; }
        }
        return b;
    }

    // -------------------- physics --------------------
    static double sellLNno(double l) {
        double L = l * l;
        return Math.sqrt(4.9048 + 0.11775 / (L - 0.0475) - 0.027169 * L);
    }
    static double sellLNne(double l) {
        double L = l * l;
        return Math.sqrt(4.5820 + 0.09921 / (L - 0.04443) - 0.021950 * L);
    }
    static double sellSiO2(double l) {
        double L = l * l;
        return Math.sqrt(1
            + 0.6961663 * L / (L - 0.0684043 * 0.0684043)
            + 0.4079426 * L / (L - 0.1162414 * 0.1162414)
            + 0.8974794 * L / (L - 9.896161 * 9.896161));
    }
    static Double slabNeff(double nC, double nS, double nCl, double d, double lam) {
        double k0 = 2 * Math.PI / lam;
        double nLo = Math.max(nS, nCl);
        if (nC <= nLo) return null;
        double dn2 = nC * nC - nS * nS;
        if (dn2 <= 0) return null;
        double V = k0 * d * Math.sqrt(nC * nC - nLo * nLo);
        if (V < 0.01) return null;
        final double a = (nS * nS - nCl * nCl) / dn2;
        final double Vs = k0 * d * Math.sqrt(dn2);
        Fn eig = b -> {
            if (b <= 0 || b >= 1) return 1e6;
            double sq = Math.sqrt(1 - b);
            double rhs = Math.atan(Math.sqrt(b / (1 - b)));
            double a2 = (b + a) / (1 - b);
            if (a2 < 0) return 1e6;
            rhs += Math.atan(Math.sqrt(a2));
            return Vs * sq - rhs;
        };
        if (eig.f(1e-12) * eig.f(1 - 1e-12) > 0) return nS;
        Double b = brentq(eig, 1e-12, 1 - 1e-12, 1e-14, 200);
        if (b == null) return nS;
        return Math.sqrt(nS * nS + b * dn2);
    }
    /** @return {Z0, C, L} for the quasi-static CPW. */
    static double[] cpwQuasiStatic(double w, double s, double eff) {
        double k = w / (w + 2 * s);
        double kp = Math.sqrt(1 - k * k);
        double Z0 = (30 * Math.PI / Math.sqrt(eff)) * (ellipk(kp * kp) / ellipk(k * k));
        double vp = C0 / Math.sqrt(eff);
        return new double[] { Z0, 1 / (Z0 * vp), Z0 / vp };
    }
    static double effEps(double tox, double s, double eS, double eL) {
        double eta = 1 - Math.exp(-2 * tox / s);
        return 0.5 * (1 + eta * eS + (1 - eta) * eL);
    }
    static double Rpm(double f, double w, double t, double sig) {
        return 2.2 * Math.sqrt(Math.PI * f * MU0 / sig) / (w + 2 * t);
    }
    static double Gpm(double f, double C, double tS, double tL, double tox, double s) {
        double eta = 1 - Math.exp(-2 * tox / s);
        return 2 * Math.PI * f * C * (eta * tS + (1 - eta) * tL);
    }
    /** @return {S11, S21, Zc, gamma}. */
    static Cx[] tlineS(double R, double L, double G, double C, double f, double len, double zref) {
        double w = 2 * Math.PI * f;
        Cx num = new Cx(R, w * L), den = new Cx(G, w * C);
        Cx gamma = num.mul(den).sqrt();
        Cx zc = num.div(den).sqrt();
        Cx gl = gamma.mul(new Cx(len, 0));
        Cx A = gl.cosh();
        Cx B = zc.mul(gl.sinh());
        Cx Cm = new Cx(1, 0).div(zc).mul(gl.sinh());
        Cx D = gl.cosh();
        Cx zr = new Cx(zref, 0);
        Cx den2 = A.add(B.div(zr)).add(Cm.mul(zr)).add(D);
        Cx s11 = A.add(B.div(zr)).sub(Cm.mul(zr)).sub(D).div(den2);
        Cx s21 = new Cx(2, 0).div(den2);
        return new Cx[] { s11, s21, zc, gamma };
    }
    static double vpiDC(double lam, double gap, double n, double r, double G, double Lel, boolean pp) {
        return lam * gap / ((pp ? 2 : 1) * n * n * n * r * G * Lel);
    }
    static double twEO(double f, double alpha, double nRF, double nOpt, double Lel) {
        double dBeta = 2 * Math.PI * f / C0 * (nRF - nOpt);
        Cx uL = new Cx(alpha, dBeta).mul(new Cx(Lel, 0));
        if (uL.abs() < 1e-12) return 1;
        Cx one = new Cx(1, 0);
        return one.sub(new Cx(-uL.re, -uL.im).exp()).div(uL).abs();
    }

    // -------------------- driver --------------------
    static String fg(double v) {
        return Double.isNaN(v) ? "N/A" : String.format("%.2f GHz", v / 1e9);
    }

    public static void main(String[] args) {
        Map<String, String> kv = new HashMap<>();
        for (String a : args) {
            int eq = a.indexOf('=');
            if (eq > 0) kv.put(a.substring(0, eq), a.substring(eq + 1));
        }
        java.util.function.BiFunction<String, Double, Double> d =
            (k, def) -> kv.containsKey(k) ? Double.parseDouble(kv.get(k)) : def;

        double w = d.apply("w", 8.0) * 1e-6;
        double s = d.apply("s", 16.0) * 1e-6;
        double tAu = d.apply("t_au", 20.0) * 1e-6;
        double tox = d.apply("tox", 1.5) * 1e-6;
        double lineLen = d.apply("line_len", 20.0) * 1e-3;
        double sigmaAu = d.apply("sigma_au", 4.1e7);
        double epsrSiO2 = d.apply("epsr_sio2", 3.9);
        double epsrLN = d.apply("epsr_ln", 28.0);
        double tanSiO2 = d.apply("tan_sio2", 0.0005);
        double tanLN = d.apply("tan_ln", 0.003);
        double lamUm = d.apply("lam", 1.55);
        double lam = lamUm * 1e-6;
        double dnTi = d.apply("dn_ti", 0.008);
        double dTi = d.apply("d_ti", 4.0) * 1e-6;
        String pol = kv.getOrDefault("pol", "TM");
        double r33 = d.apply("r33", 30.8) * 1e-12;
        double r13 = d.apply("r13", 8.6) * 1e-12;
        double gammaOv = d.apply("gamma", 0.5);
        boolean pp = !"0".equals(kv.getOrDefault("pp", "1"));
        double fStart = d.apply("f_start", 1.0) * 1e9;
        double fStop = d.apply("f_stop", 67.0) * 1e9;
        int npts = Math.max(2, (int) (double) d.apply("npts", 800.0));
        double z0Port = d.apply("z0_port", 50.0);

        double[] f = linspace(fStart, fStop, npts);

        double epsRf = effEps(tox, s, epsrSiO2, epsrLN);
        double[] qs = cpwQuasiStatic(w, s, epsRf);
        double Z0q = qs[0], Cpm = qs[1], Lpm = qs[2];

        double[] s21dB = new double[f.length], alpha = new double[f.length];
        double[] zcRe = new double[f.length], nRf = new double[f.length];
        double[] gReal = new double[f.length];
        for (int i = 0; i < f.length; i++) {
            double R = Rpm(f[i], w, tAu, sigmaAu);
            double G = Gpm(f[i], Cpm, tanSiO2, tanLN, tox, s);
            Cx[] tl = tlineS(R, Lpm, G, Cpm, f[i], lineLen, z0Port);
            s21dB[i] = 20 * Math.log10(Math.max(tl[1].abs(), 1e-15));
            alpha[i] = tl[3].re * 8.686 / 100;
            zcRe[i] = tl[2].re;
            gReal[i] = tl[3].re;
            nRf[i] = C0 * tl[3].im / (2 * Math.PI * f[i]);
        }
        double f3 = findCrossing(f, s21dB, -3);
        double f6 = findCrossing(f, s21dB, -6);
        double nRfQs = Math.sqrt(epsRf);

        double noB = sellLNno(lamUm), neB = sellLNne(lamUm), nSiO2 = sellSiO2(lamUm);
        boolean isTM = pol.equals("TM");
        double nBulk = isTM ? neB : noB;
        String polTag = isTM ? "TM (ne)" : "TE (no)";
        double nCore = nBulk + dnTi;
        Double nOptB = slabNeff(nCore, nBulk, nSiO2, dTi, lam);
        double nOpt = nOptB == null ? nBulk : nOptB;
        double dl = 0.001;
        double dnDl = ((isTM ? sellLNne(lamUm + dl) : sellLNno(lamUm + dl))
                     - (isTM ? sellLNne(lamUm - dl) : sellLNno(lamUm - dl))) / (2 * dl);
        double ng = nBulk - lamUm * dnDl;

        double rEO = isTM ? r33 : r13;
        double nEO = isTM ? neB : noB;
        double gap = s;
        double vpiDc = vpiDC(lam, gap, nEO, rEO, gammaOv, lineLen, pp);
        double vpilDc = vpiDc * lineLen * 100;

        double[] mEoDB = new double[f.length], vpif = new double[f.length];
        for (int i = 0; i < f.length; i++) {
            double m = Math.max(twEO(f[i], gReal[i], nRf[i], nOpt, lineLen), 1e-15);
            mEoDB[i] = 20 * Math.log10(m);
            vpif[i] = vpiDc / m;
        }
        double fEo3 = findCrossing(f, mEoDB, -3);
        double fEo6 = findCrossing(f, mEoDB, -6);

        StringBuilder o = new StringBuilder();
        o.append("== RF (CPW) ==========================\n");
        o.append(String.format("  eps_eff (RF)   = %.4f%n", epsRf));
        o.append(String.format("  n_RF (q.s.)    = %.4f%n", nRfQs));
        o.append(String.format("  Z0 (q.s.)      = %.2f ohm%n", Z0q));
        o.append(String.format("  L              = %.4f uH/m%n", Lpm * 1e6));
        o.append(String.format("  C              = %.4f pF/m%n", Cpm * 1e12));
        o.append(String.format("  v_ph (RF)      = %.4f x10^8 m/s%n", C0 / nRfQs / 1e8));
        o.append("\n");
        o.append(String.format("  -3 dB freq     = %s%n", fg(f3)));
        o.append(String.format("  -6 dB freq     = %s%n", fg(f6)));
        o.append(String.format("  alpha @ 10 GHz = %.4f dB/cm%n", interp(10e9, f, alpha)));
        o.append(String.format("  alpha @ 40 GHz = %.4f dB/cm%n", interp(40e9, f, alpha)));
        o.append(String.format("  Re(Zc)@10G     = %.2f ohm%n", interp(10e9, f, zcRe)));
        o.append("\n");
        o.append("== Optical (LiNbO3) ==================\n");
        o.append(String.format("  lambda         = %s um%n", lamUm));
        o.append(String.format("  Pol.           = %s%n", polTag));
        o.append(String.format("  no (bulk)      = %.5f%n", noB));
        o.append(String.format("  ne (bulk)      = %.5f%n", neB));
        o.append(String.format("  n_SiO2         = %.5f%n", nSiO2));
        o.append(String.format("  n_core (+dn)   = %.5f%n", nCore));
        o.append(String.format("  n_eff (WG)     = %.5f%n", nOpt));
        o.append(String.format("  n_g (group)    = %.5f%n", ng));
        o.append("\n");
        o.append("== Velocity Matching =================\n");
        double dnMid = interp(0.5 * (fStart + fStop), f, nRf) - nOpt;
        o.append(String.format("  dn (mid-freq)  = %+.5f%n", dnMid));
        o.append(String.format("  n_RF - n_opt   = %+.5f  (q.s.)%n", nRfQs - nOpt));
        o.append("\n");
        o.append("== Vpi / EO ==========================\n");
        o.append(String.format("  EO coeff       = %s = %.1f pm/V%n", isTM ? "r33" : "r13", rEO * 1e12));
        o.append(String.format("  n (EO)         = %.5f%n", nEO));
        o.append(String.format("  Gamma (overlap)= %.3f%n", gammaOv));
        o.append(String.format("  MZM type       = %s%n", pp ? "Push-pull" : "Single-arm"));
        o.append(String.format("  Electrode gap  = %.1f um%n", gap * 1e6));
        o.append(String.format("  Electrode L    = %.1f mm%n", lineLen * 1e3));
        o.append("\n");
        o.append(String.format("  Vpi (DC)       = %.3f V%n", vpiDc));
        o.append(String.format("  Vpi.L (DC)     = %.3f V.cm%n", vpilDc));
        o.append("\n");
        o.append(String.format("  Vpi @ 10 GHz   = %.3f V%n", interp(10e9, f, vpif)));
        o.append(String.format("  Vpi @ 40 GHz   = %.3f V%n", interp(40e9, f, vpif)));
        o.append("\n");
        o.append(String.format("  EO BW (-3dB el)= %s%n", fg(fEo3)));
        o.append(String.format("  EO BW (-6dB op)= %s%n", fg(fEo6)));
        System.out.print(o);
    }
}
