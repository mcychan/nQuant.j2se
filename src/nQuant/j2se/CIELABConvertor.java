package nQuant.j2se;

import java.awt.Color;
import java.math.BigDecimal;

public class CIELABConvertor {
	private final static char BYTE_MAX = -Byte.MIN_VALUE + Byte.MAX_VALUE;
	private final static double M_PI = 3.14159265358979323846;
	
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
	
	static Lab RGB2LAB(final int pixel)
	{
		final Color c1 = new Color(pixel);
		float r = c1.getRed() / 255f, g = c1.getGreen() / 255f, b = c1.getBlue() / 255f;
		float x, y, z;

		r = (r > 0.04045) ? (float) Math.pow((r + 0.055) / 1.055, 2.4) : r / 12.92f;
		g = (g > 0.04045) ? (float) Math.pow((g + 0.055) / 1.055, 2.4) : g / 12.92f;
		b = (b > 0.04045) ? (float) Math.pow((b + 0.055) / 1.055, 2.4) : b / 12.92f;

		x = (r * 0.4124f + g * 0.3576f + b * 0.1805f) / 0.95047f;
		y = (r * 0.2126f + g * 0.7152f + b * 0.0722f) / 1.00000f;
		z = (r * 0.0193f + g * 0.1192f + b * 0.9505f) / 1.08883f;

		x = (x > 0.008856) ? (float) Math.cbrt(x) : (7.787f * x) + 16f / 116f;
		y = (y > 0.008856) ? (float) Math.cbrt(y) : (7.787f * y) + 16f / 116f;
		z = (z > 0.008856) ? (float) Math.cbrt(z) : (7.787f * z) + 16f / 116f;

		Lab lab = new Lab();
		lab.alpha = c1.getAlpha();
		lab.L = (116f * y) - 16;
		lab.A = 500f * (x - y);
		lab.B = 200f * (y - z);
		return lab;
	}
		
	static Color LAB2RGB(final Lab lab){
		float y = (lab.L + 16) / 116;
		float x = lab.A / 500 + y;
		float z = y - lab.B / 200;
		float r, g, b;

		x = 0.95047f * ((x * x * x > 0.008856) ? x * x * x : (x - 16f / 116f) / 7.787f);
		y = 1.00000f * ((y * y * y > 0.008856) ? y * y * y : (y - 16f / 116f) / 7.787f);
		z = 1.08883f * ((z * z * z > 0.008856) ? z * z * z : (z - 16f / 116f) / 7.787f);

		r = x *  3.2406f + y * -1.5372f + z * -0.4986f;
		g = x * -0.9689f + y *  1.8758f + z *  0.0415f;
		b = x *  0.0557f + y * -0.2040f + z *  1.0570f;

		r = (r > 0.0031308) ? (float) (1.055f * Math.pow(r, 1.0 / 2.4) - 0.055) : 12.92f * r;
		g = (g > 0.0031308) ? (float) (1.055f * Math.pow(g, 1.0 / 2.4) - 0.055) : 12.92f * g;
		b = (b > 0.0031308) ? (float) (1.055f * Math.pow(b, 1.0 / 2.4) - 0.055) : 12.92f * b;

		return new Color((float) Math.max(0, Math.min(1, r)), (float) Math.max(0, Math.min(1, g)), (float) Math.max(0, Math.min(1, b)), (float) Math.min(lab.alpha, BYTE_MAX) / 255.0f);
	}

	/*******************************************************************************
	* Conversions.
	******************************************************************************/

	private static final float deg2Rad(final float deg)
	{
		return (float) (deg * (M_PI / 180f));
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
}
