/**
 * Refractive index of silicon nitride (Si3N4) via the Sellmeier equation.
 *
 * Luke et al., Opt. Lett. 40, 4823 (2015) — LPCVD stoichiometric Si3N4:
 *
 *   n(lambda)^2 = 1 + B1*lambda^2/(lambda^2 - C1) + B2*lambda^2/(lambda^2 - C2)
 *
 *     B1 = 3.0249,   C1 = 0.1353406^2
 *     B2 = 40314.0,  C2 = 1239.842^2
 *
 * Also reports the group index n_g = n - lambda * dn/dlambda.
 *
 * Valid range: 0.31 - 5.504 um.
 *
 * Compile: javac NSi3N4.java
 * Run:     java NSi3N4 1.55           (one or more wavelengths in um)
 *          java NSi3N4 1.31 1.55 0.85
 */
public class NSi3N4 {

    static final double B1 = 3.0249, B2 = 40314.0;
    static final double C1 = 0.1353406 * 0.1353406;
    static final double C2 = 1239.842 * 1239.842;

    public static double nSi3N4(double wavelengthUm) {
        double l2 = wavelengthUm * wavelengthUm;
        double n2 = 1.0 + B1 * l2 / (l2 - C1) + B2 * l2 / (l2 - C2);
        return Math.sqrt(n2);
    }

    /**
     * Group index n_g = n - lambda * dn/dlambda. From the Sellmeier equation:
     *   n_g = n + (lambda^2 / n) * sum_i B_i * C_i / (lambda^2 - C_i)^2.
     */
    public static double nGroup(double wavelengthUm) {
        double l2 = wavelengthUm * wavelengthUm;
        double n = nSi3N4(wavelengthUm);
        double s = B1 * C1 / Math.pow(l2 - C1, 2)
                 + B2 * C2 / Math.pow(l2 - C2, 2);
        return n + (l2 / n) * s;
    }

    public static void main(String[] args) {
        if (args.length == 0) {
            System.out.println("Usage: java NSi3N4 <wavelength_um> [<wavelength_um> ...]");
            System.out.println("Example: java NSi3N4 1.55");
            return;
        }
        for (String arg : args) {
            double wl = Double.parseDouble(arg);
            System.out.printf("lambda = %.4f um   n = %.6f   n_g = %.6f%n", wl, nSi3N4(wl), nGroup(wl));
        }
    }
}
