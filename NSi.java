/**
 * Refractive index of crystalline silicon (Si) via the Sellmeier equation.
 *
 * Salzberg & Villa (1957), as tabulated for crystalline silicon:
 *
 *   n(lambda)^2 = 1 + B1*l^2/(l^2 - C1) + B2*l^2/(l^2 - C2) + B3*l^2/(l^2 - C3)
 *
 *     B1 = 10.6684293,    C1 = 0.301516485^2
 *     B2 = 0.0030434748,  C2 = 1.13475115^2
 *     B3 = 1.54133408,    C3 = 1104.0^2
 *
 * Also reports the group index n_g = n - lambda * dn/dlambda.
 *
 * Valid range: 1.36 - 11 um.
 *
 * Compile: javac NSi.java
 * Run:     java NSi 1.55           (one or more wavelengths in um)
 *          java NSi 1.31 1.55 3.0
 */
public class NSi {

    static final double B1 = 10.6684293, B2 = 0.0030434748, B3 = 1.54133408;
    static final double C1 = 0.301516485 * 0.301516485;
    static final double C2 = 1.13475115 * 1.13475115;
    static final double C3 = 1104.0 * 1104.0;

    public static double nSi(double wavelengthUm) {
        double l2 = wavelengthUm * wavelengthUm;
        double n2 = 1.0 + B1 * l2 / (l2 - C1) + B2 * l2 / (l2 - C2) + B3 * l2 / (l2 - C3);
        return Math.sqrt(n2);
    }

    /**
     * Group index n_g = n - lambda * dn/dlambda. From the Sellmeier equation:
     *   n_g = n + (lambda^2 / n) * sum_i B_i * C_i / (lambda^2 - C_i)^2.
     */
    public static double nGroup(double wavelengthUm) {
        double l2 = wavelengthUm * wavelengthUm;
        double n = nSi(wavelengthUm);
        double s = B1 * C1 / Math.pow(l2 - C1, 2)
                 + B2 * C2 / Math.pow(l2 - C2, 2)
                 + B3 * C3 / Math.pow(l2 - C3, 2);
        return n + (l2 / n) * s;
    }

    public static void main(String[] args) {
        if (args.length == 0) {
            System.out.println("Usage: java NSi <wavelength_um> [<wavelength_um> ...]");
            System.out.println("Example: java NSi 1.55");
            return;
        }
        for (String arg : args) {
            double wl = Double.parseDouble(arg);
            System.out.printf("lambda = %.4f um   n = %.6f   n_g = %.6f%n", wl, nSi(wl), nGroup(wl));
        }
    }
}
