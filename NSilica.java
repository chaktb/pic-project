/**
 * Refractive index of fused silica via the Sellmeier equation (Malitson 1965).
 *
 *   n(lambda)^2 = 1 + sum_i  B_i * lambda^2 / (lambda^2 - C_i)
 *
 *     B1 = 0.6961663,  C1 = 0.0684043^2
 *     B2 = 0.4079426,  C2 = 0.1162414^2
 *     B3 = 0.8974794,  C3 = 9.896161^2
 *
 * Valid range: 0.21 - 6.7 um.
 *
 * Compile: javac NSilica.java
 * Run:     java NSilica 1.55           (one or more wavelengths in um)
 *          java NSilica 1.31 1.55 0.85
 */
public class NSilica {

    public static double nSilica(double wavelengthUm) {
        double B1 = 0.6961663, B2 = 0.4079426, B3 = 0.8974794;
        double C1 = 0.0684043 * 0.0684043;
        double C2 = 0.1162414 * 0.1162414;
        double C3 = 9.896161 * 9.896161;
        double l2 = wavelengthUm * wavelengthUm;
        double n2 = 1.0 + B1 * l2 / (l2 - C1) + B2 * l2 / (l2 - C2) + B3 * l2 / (l2 - C3);
        return Math.sqrt(n2);
    }

    public static void main(String[] args) {
        if (args.length == 0) {
            System.out.println("Usage: java NSilica <wavelength_um> [<wavelength_um> ...]");
            System.out.println("Example: java NSilica 1.55");
            return;
        }
        for (String arg : args) {
            double wl = Double.parseDouble(arg);
            System.out.printf("lambda = %.4f um   n = %.6f%n", wl, nSilica(wl));
        }
    }
}
