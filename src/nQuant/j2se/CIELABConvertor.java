package nQuant.j2se;

import java.awt.Color;
import java.math.BigDecimal;

public class CIELABConvertor {
	protected static final char BYTE_MAX = -Byte.MIN_VALUE + Byte.MAX_VALUE;
	private static final double XYZ_WHITE_REFERENCE_X = 95.047;
	private static final double XYZ_WHITE_REFERENCE_Y = 100;
	private static final double XYZ_WHITE_REFERENCE_Z = 108.883;
	private static final double XYZ_EPSILON = 0.008856;
	private static final double XYZ_KAPPA = 903.3;
	
	static class MutableDouble extends Number {

		private static final long serialVersionUID = -8826262264116498065L;
		private double value;

		public MutableDouble(double value) {
			this.value = value;
		}
		
		public MutableDouble() {
			this(0.0);
		}		

		public void setValue(double value) {
			this.value = value;
		}

		@Override
		public int intValue() {
			return (int) value;
		}

		@Override
		public long longValue() {
			return (long) value;
		}

		@Override
		public float floatValue() {
			return (float) value;
		}

		@Override
		public double doubleValue() {
			return value;
		}
		
	}
	
	static class Lab {
		float alpha = BYTE_MAX;
		float A = 0.0f;
		float B = 0.0f;
		float L = 0.0f;
	}
	
	static int clamp(double val, int min, int max) {
		return (int) Math.max(min, Math.min(max, val));
	}
	
	private static float pivotXyzComponent(double component) {
		return component > XYZ_EPSILON
				? (float) Math.cbrt(component)
				: (float) ((XYZ_KAPPA * component + 16) / 116.0);
	}
	
	protected static double gammaToLinear(int channel)
	{
		final double c = channel / 255.0;
		return c < 0.04045 ? c / 12.92 : Math.pow((c + 0.055) / 1.055, 2.4);
	}
	
	static Lab RGB2LAB(final int pixel, final boolean hasAlpha)
	{
		final Color c1 = new Color(pixel, hasAlpha);
		double sr = gammaToLinear(c1.getRed());
		double sg = gammaToLinear(c1.getGreen());
		double sb = gammaToLinear(c1.getBlue());
		final float x = pivotXyzComponent(100 * (sr * 0.4124 + sg * 0.3576 + sb * 0.1805) / XYZ_WHITE_REFERENCE_X);
		final float y = pivotXyzComponent(100 * (sr * 0.2126 + sg * 0.7152 + sb * 0.0722) / XYZ_WHITE_REFERENCE_Y);
		final float z = pivotXyzComponent(100 * (sr * 0.0193 + sg * 0.1192 + sb * 0.9505) / XYZ_WHITE_REFERENCE_Z);        

		Lab lab = new Lab();
		lab.alpha = c1.getAlpha();
		lab.L = Math.max(0, 116 * y - 16);
		lab.A = 500 * (x - y);
		lab.B = 200 * (y - z);
		return lab;
	}
		
	static Color LAB2RGB(final Lab lab)
	{
		final double fy = (lab.L + 16.0) / 116.0;
		final double fx = lab.A / 500 + fy;
		final double fz = fy - lab.B / 200.0;
		double tmp = fx * fx * fx;
		final double xr = tmp > XYZ_EPSILON ? tmp : (116.0 * fx - 16) / XYZ_KAPPA;
		final double yr = lab.L > XYZ_KAPPA * XYZ_EPSILON ? fy * fy * fy : lab.L / XYZ_KAPPA;
		tmp = fz * fz * fz;
		final double zr = tmp > XYZ_EPSILON ? tmp : (116.0 * fz - 16) / XYZ_KAPPA;
		final double x = xr * XYZ_WHITE_REFERENCE_X;
		final double y = yr * XYZ_WHITE_REFERENCE_Y;
		final double z = zr * XYZ_WHITE_REFERENCE_Z;
		
		int a = clamp(lab.alpha, 0, 255);
		double r = (x * 3.2406 + y * -1.5372 + z * -0.4986) / 100.0;
		double g = (x * -0.9689 + y * 1.8758 + z * 0.0415) / 100.0;
		double b = (x * 0.0557 + y * -0.2040 + z * 1.0570) / 100.0;
		r = r > 0.0031308 ? 1.055 * Math.pow(r, 1 / 2.4) - 0.055 : 12.92 * r;
		g = g > 0.0031308 ? 1.055 * Math.pow(g, 1 / 2.4) - 0.055 : 12.92 * g;
		b = b > 0.0031308 ? 1.055 * Math.pow(b, 1 / 2.4) - 0.055 : 12.92 * b;

		return new Color((int) clamp(Math.rint(r * BYTE_MAX), 0, BYTE_MAX),
			(int) clamp(Math.rint(g * BYTE_MAX), 0, BYTE_MAX),
			(int) clamp(Math.rint(b * BYTE_MAX), 0, BYTE_MAX), a);
	}

	/*******************************************************************************
	* Conversions.
	******************************************************************************/

	private static final float deg2Rad(final float deg)
	{
		return (float) (deg * (Math.PI / 180f));
	}

	static float L_prime_div_k_L_S_L(final Lab lab1, final Lab lab2)
	{
		final float k_L = 1.0f;
		float deltaLPrime = lab2.L - lab1.L;
		float barLPrime = (lab1.L + lab2.L) / 2f;
		float S_L = (float)(1 + ((0.015f * Math.pow(barLPrime - 50f, 2f)) / Math.sqrt(20 + Math.pow(barLPrime - 50f, 2f))));
		return deltaLPrime / (k_L * S_L);
	}

	static float C_prime_div_k_L_S_L(final Lab lab1, final Lab lab2, MutableDouble a1Prime, MutableDouble a2Prime, MutableDouble CPrime1, MutableDouble CPrime2)
	{
		final float k_C = 1f;
		final float pow25To7 = 6103515625f; /* pow(25, 7) */
		float C1 = (float)(Math.sqrt((lab1.A * lab1.A) + (lab1.B * lab1.B)));
		float C2 = (float)(Math.sqrt((lab2.A * lab2.A) + (lab2.B * lab2.B)));
		float barC = (C1 + C2) / 2f;
		float G = (float)(0.5f * (1 - Math.sqrt(Math.pow(barC, 7) / (Math.pow(barC, 7) + pow25To7))));
		a1Prime.setValue((1.0 + G) * lab1.A);
		a2Prime.setValue((1.0 + G) * lab2.A);

		CPrime1.setValue(Math.sqrt((a1Prime.doubleValue() * a1Prime.doubleValue()) + (lab1.B * lab1.B)));
		CPrime2.setValue(Math.sqrt((a2Prime.doubleValue() * a2Prime.doubleValue()) + (lab2.B * lab2.B)));
		float deltaCPrime = CPrime2.floatValue() - CPrime1.floatValue();
		float barCPrime = (CPrime1.floatValue() + CPrime2.floatValue()) / 2f;
		
		float S_C = 1 + (0.045f * barCPrime);
		return deltaCPrime / (k_C * S_C);
	}

	static float H_prime_div_k_L_S_L(final Lab lab1, final Lab lab2, final Number a1Prime, final Number a2Prime, final Number CPrime1, final Number CPrime2, MutableDouble barCPrime, MutableDouble barhPrime)
	{
		final float k_H = 1f;
		final float deg360InRad = deg2Rad(360f);
		final float deg180InRad = deg2Rad(180f);
		double CPrimeProduct = CPrime1.doubleValue() * CPrime2.doubleValue();
		double hPrime1;
		if (BigDecimal.ZERO.equals(new BigDecimal(lab1.B)) && BigDecimal.ZERO.equals(new BigDecimal(a1Prime.doubleValue())))
			hPrime1 = 0.0;
		else {
			hPrime1 = Math.atan2(lab1.B, a1Prime.doubleValue());
			/*
			* This must be converted to a hue angle in degrees between 0
			* and 360 by addition of 2π to negative hue angles.
			*/
			if (hPrime1 < 0)
				hPrime1 += deg360InRad;
		}
		double hPrime2;
		if (BigDecimal.ZERO.equals(new BigDecimal(lab2.B)) && BigDecimal.ZERO.equals(new BigDecimal(a2Prime.doubleValue())))
			hPrime2 = 0.0;
		else {
			hPrime2 = Math.atan2(lab2.B, a2Prime.doubleValue());
			/*
			* This must be converted to a hue angle in degrees between 0
			* and 360 by addition of 2π to negative hue angles.
			*/
			if (hPrime2 < 0)
				hPrime2 += deg360InRad;
		}
		double deltahPrime;
		if (BigDecimal.ZERO.equals(new BigDecimal(CPrimeProduct)))
			deltahPrime = 0;
		else {
			/* Avoid the Math.abs() call */
			deltahPrime = hPrime2 - hPrime1;
			if (deltahPrime < -deg180InRad)
				deltahPrime += deg360InRad;
			else if (deltahPrime > deg180InRad)
				deltahPrime -= deg360InRad;
		}

		double deltaHPrime = 2.0 * Math.sqrt(CPrimeProduct) * Math.sin(deltahPrime / 2.0);
		double hPrimeSum = hPrime1 + hPrime2;
		if (BigDecimal.ZERO.equals(new BigDecimal(CPrime1.doubleValue() * CPrime2.doubleValue()))) {
			barhPrime.setValue(hPrimeSum);
		}
		else {
			if (Math.abs(hPrime1 - hPrime2) <= deg180InRad)
				barhPrime.setValue(hPrimeSum / 2.0);
			else {
				if (hPrimeSum < deg360InRad)
					barhPrime.setValue((hPrimeSum + deg360InRad) / 2.0);
				else
					barhPrime.setValue((hPrimeSum - deg360InRad) / 2.0);
			}
		}

		barCPrime.setValue((CPrime1.doubleValue() + CPrime2.doubleValue()) / 2.0);
		double T = 1.0 - (0.17 * Math.cos(barhPrime.doubleValue() - deg2Rad(30f))) +
			(0.24 * Math.cos(2.0 * barhPrime.doubleValue())) +
			(0.32 * Math.cos((3.0 * barhPrime.doubleValue()) + deg2Rad(6f))) -
			(0.20 * Math.cos((4.0 * barhPrime.doubleValue()) - deg2Rad(63f)));
		double S_H = 1 + (0.015f * barCPrime.doubleValue() * T);
		return (float) (deltaHPrime / (k_H * S_H));
	}

	static float R_T(final Number barCPrime, final Number barhPrime, final float C_prime_div_k_L_S_L, final float H_prime_div_k_L_S_L)
	{
		final double pow25To7 = 6103515625.0; /* Math.pow(25, 7) */
		double deltaTheta = deg2Rad(30f) * Math.exp(-Math.pow((barhPrime.doubleValue() - deg2Rad(275f)) / deg2Rad(25f), 2.0));
		double R_C = 2.0 * Math.sqrt(Math.pow(barCPrime.doubleValue(), 7.0) / (Math.pow(barCPrime.doubleValue(), 7.0) + pow25To7));
		double R_T = (-Math.sin(2.0 * deltaTheta)) * R_C;
		return (float) (R_T * C_prime_div_k_L_S_L * H_prime_div_k_L_S_L);
	}

	/* From the paper "The CIEDE2000 Color-Difference Formula: Implementation Notes, */
	/* Supplementary Test Data, and Mathematical Observations", by */
	/* Gaurav Sharma, Wencheng Wu and Edul N. Dalal, */
	/* Color Res. Appl., vol. 30, no. 1, pp. 21-30, Feb. 2005. */
	/* Return the CIEDE2000 Delta E color difference measure squared, for two Lab values */
	static float CIEDE2000(final Lab lab1, final Lab lab2)
	{
		float deltaL_prime_div_k_L_S_L = L_prime_div_k_L_S_L(lab1, lab2);
		MutableDouble a1Prime = new MutableDouble(), a2Prime = new MutableDouble(), CPrime1 = new MutableDouble(), CPrime2 = new MutableDouble();
		float deltaC_prime_div_k_L_S_L = C_prime_div_k_L_S_L(lab1, lab2, a1Prime, a2Prime, CPrime1, CPrime2);
		MutableDouble barCPrime = new MutableDouble(), barhPrime = new MutableDouble();
		float deltaH_prime_div_k_L_S_L = H_prime_div_k_L_S_L(lab1, lab2, a1Prime, a2Prime, CPrime1, CPrime2, barCPrime, barhPrime);
		float deltaR_T = R_T(barCPrime, barhPrime, deltaC_prime_div_k_L_S_L, deltaH_prime_div_k_L_S_L);
		return (float) (Math.pow(deltaL_prime_div_k_L_S_L, 2.0) +
			Math.pow(deltaC_prime_div_k_L_S_L, 2.0) +
			Math.pow(deltaH_prime_div_k_L_S_L, 2.0) +
			deltaR_T);
	}
	
	static double Y_Diff(final Color c1, final Color c2)
	{
		java.util.function.Function<Color, Double> color2Y = c -> {
			double sr = gammaToLinear(c.getRed());
			double sg = gammaToLinear(c.getGreen());
			double sb = gammaToLinear(c.getBlue());
			return sr * 0.2126 + sg * 0.7152 + sb * 0.0722;
		};
		
		double y = color2Y.apply(c1);
		double y2 = color2Y.apply(c2);
		return Math.abs(y2 - y) * XYZ_WHITE_REFERENCE_Y;
	}
}
