package nQuant.j2se;

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.awt.image.ColorModel;
import java.awt.image.DataBuffer;
import java.awt.image.DirectColorModel;
import java.awt.image.ImageObserver;
import java.awt.image.IndexColorModel;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;

import nQuant.j2se.CIELABConvertor.Lab;
import nQuant.j2se.CIELABConvertor.MutableDouble;

import static nQuant.j2se.BitmapUtilities.BYTE_MAX;

public class NeuQuantizer {
	//====================
	// NeuralNet Color quantizing

	/* NeuQuant Neural-Net Quantization Algorithm
	* ------------------------------------------
	*
	* Copyright (c) 1994 Anthony Dekker
	* Copyright (c) 2019-2021 Miller Cy Chan
	*
	* NEUQUANT Neural-Net quantization algorithm by Anthony Dekker, 1994.
	* See "Kohonen neural networks for optimal colour quantization"
	* in "Network: Computation in Neural Systems" Vol. 5 (1994) pp 351-367.
	* for research paper of the algorithm.
	* See also  https://www.researchgate.net/publication/232079905_Kohonen_neural_networks_for_optimal_colour_quantization
	*
	* Any party obtaining a copy of these files from the author, directly or
	* indirectly, is granted, free of charge, a full and unrestricted irrevocable,
	* world-wide, paid up, royalty-free, nonexclusive right and license to deal
	* in this software and documentation files (the "Software"), including without
	* limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
	* and/or sell copies of the Software, and to permit persons who receive
	* copies from any such party to do so, with the only requirement being
	* that this copyright notice remain intact.
	*/

	private byte alphaThreshold = 0xF;
	private final int ncycles = 115;			// no. of learning cycles
	private final int radiusbiasshift = 8;
	private final int radiusbias = 1 << radiusbiasshift;

	private int netsize = 256;		// number of colours used	
	private int maxnetpos = netsize - 1;
	private int initrad = netsize >> 3;   // for 256 cols, radius starts at 32
	private double initradius = initrad * 1.0;
	private final int radiusdec = 30; // factor of 1/30 each cycle

	private final short normal_learning_extension_factor = 2; /* normally learn twice as long */
	private final short extra_long_colour_threshold = 40; /* learn even longer when under 40 */
	private final short extra_long_divisor = 8; /* fraction of netsize for extra learning */

	private final int alphabiasshift = 10;			// alpha starts at 1
	private final int initalpha = 1 << alphabiasshift; // biased by 10 bits
	private final int alpharadbshift = alphabiasshift + radiusbiasshift;
	private final double alpharadbias = (double)(1 << alpharadbshift);
	private final double exclusion_threshold = 0.5;

	private final short REPEL_THRESHOLD = 16;          /* See repel_coincident()... */
	private final short REPEL_STEP_DOWN = 1;              /* ... for an explanation of... */
	private final short REPEL_STEP_UP = 4;                 /* ... how these points work. */
	private int[] repel_points;

	/* defs for freq and bias */
	private final int gammashift = 10;                  /* gamma = 1024 */
	private final double gamma = (double)(1 << gammashift);
	private final int betashift = 10;
	private final double beta = (1.0 / (double)(1 << betashift));/* beta = 1/1024 */
	private final double betagamma = (double)(1 << (gammashift - betashift));

	private Lab[] network; // the network itself

	private int[] netindex; // for network lookup - really 256

	private double[] bias;  // bias and freq arrays for learning
	private double[] freq;
	private double[] radpower;

	protected boolean hasSemiTransparency = false;
	protected int m_transparentPixelIndex = -1;
	protected final int width, height;	
	protected int[] pixels;
	protected Color m_transparentColor = new Color(BYTE_MAX, BYTE_MAX, BYTE_MAX, 0);
	protected Map<Integer, Lab> pixelMap = new HashMap<>();
	protected Map<Integer, Short> nearestMap = new HashMap<>();
	private ColorModel m_colorModel;
	private static Random random = new Random();

	private static double colorimportance(double al)
	{
		double transparency = 1.0 - al / 255.0;
		return (1.0 - transparency * transparency);
	}

	/* posneg(x) returns +1 if x is positive or zero, or -1 otherwise.
	*
	* XXX posneg could be turned into a vector function in various ways.
	*/
	private static float posneg(double x) {
		if (x < 0)
			return -1f;
		return 1f;
	}
	
	private NeuQuantizer(BufferedImage im, int w, int h) {
		width = w;
		height = h;
		setPixels(im);
	}

	public NeuQuantizer(BufferedImage im, ImageObserver obs) {
		this(im, im.getWidth(obs), im.getHeight(obs));
	}

	private void setPixels(BufferedImage im) {
		pixels = im.getRGB(0, 0, width, height, null, 0, width);
	}

	private void setUpArrays() {
		network = new Lab[netsize];
		netindex = new int[Math.max(netsize, 256)];
		repel_points = new int[Math.max(netsize, 256)];
		bias = new double[netsize];
		freq = new double[netsize];
		radpower = new double[initrad];
			
		for (int i = 0; i < netsize; ++i) {
			network[i] = new Lab();
			network[i].L = network[i].A = network[i].B = i / netsize;

			/*  Sets alpha values at 0 for dark pixels. */
			if (i < 16)
				network[i].alpha = i * 16f;
			else
				network[i].alpha = BYTE_MAX;

			freq[i] = 1.0 / netsize;
		}
	}

	Lab getLab(final int pixel)
	{
		Lab lab1 = pixelMap.get(pixel);
		if (lab1 == null) {
			lab1 = CIELABConvertor.RGB2LAB(pixel, hasAlpha());
			pixelMap.put(pixel, lab1);
		}
		return lab1;
	}

	private static int round_biased(double temp)
	{
		if (temp < 0)
			return 0;
		temp = Math.floor((temp / 255.0 * 256.0));

		if (temp > BYTE_MAX)
			return BYTE_MAX;
		return (int)temp;
	}

	private void altersingle(double alpha, int i, int al, double L, double A, double B) {
		double colorimp = 1.0;//0.5;// + 0.7 * colorimportance(al);

		alpha /= initalpha;

		/* alter hit neuron */
		network[i].alpha -= alpha * (network[i].alpha - al);
		network[i].L -= colorimp * alpha * (network[i].L - L);
		network[i].A -= colorimp * alpha * (network[i].A - A);
		network[i].B -= colorimp * alpha * (network[i].B - B);
	}

	private void alterneigh(int rad, int i, int al, double L, double A, double B) {
		int lo = i - rad;
		if (lo < 0)
			lo = 0;
		int hi = i + rad;
		if (hi > maxnetpos)
			hi = maxnetpos;

		int j = i + 1;
		int k = i - 1;
		int m = 0;
		while ((j <= hi) || (k >= lo)) {			
			double learning_rate = radpower[m++] / alpharadbias;
			if (j <= hi) {
				network[j].alpha -= learning_rate * (network[j].alpha - al);
				network[j].L -= learning_rate * (network[j].L - L);
				network[j].A -= learning_rate * (network[j].A - A);
				network[j].B -= learning_rate * (network[j].B - B);
				j++;
			}
			if (k >= lo) {
				network[k].alpha -= learning_rate * (network[k].alpha - al);
				network[k].L -= learning_rate * (network[k].L - L);
				network[k].A -= learning_rate * (network[k].A - A);
				network[k].B -= learning_rate * (network[k].B - B);
				k--;
			}
		}
	}

	/* repelcoincident(i) traverses the entire neural network, identifies any neurons that are close in colour to neuron i, and
	 * moves them away from neuron i.  Our definition of 'close' is being within exclusion_threshold of neuron i in each component.
	 *
	 * Because repelcoincident() can be cpu-intensive, each neuron is given repel points (in the repel_points array) to limit how
	 * often a full repel pass is made.  Each time neuron i gets a full repel pass, it gets REPEL_STEP_UP points.  When the number
	 * of points it has exceeds REPEL_THRESHOLD, we stop doing full repel passes, and deduct REPEL_STEP_DOWN points instead.
	 * Eventually the number of repel points will eventually oscillate around the threshold.  With current settings, that means
	 * that only every 4th function call will result in a full pass.
 	*/
	private void repelcoincident(int i) {
		/* Use brute force to precompute the distance vectors between our neuron and each neuron. */

		if (repel_points[i] > REPEL_THRESHOLD) {
			repel_points[i] -= REPEL_STEP_DOWN;
			return;
		}

		Lab[] diffs = new Lab[netsize];
		for (int vdx = 0; vdx < netsize; ++vdx) {
			diffs[vdx] = new Lab();
			diffs[vdx].alpha = network[vdx].alpha - network[i].alpha;
			diffs[vdx].L = Math.abs(network[vdx].L - network[i].L);
			diffs[vdx].A = Math.abs(network[vdx].A - network[i].A);
			diffs[vdx].B = Math.abs(network[vdx].B - network[i].B);
		}


		/* Identify which neurons are too close to neuron[i] and shift them away.
		 * */

		 /* repel_step is the amount we move the neurons away by in each component.
		  * radpower[0]/alpharadbias is similar to the proportion alterneigh uses.
		  * */
		double repel_step = exclusion_threshold * (radpower[0] / alpharadbias);

		for (int j = 0; j < netsize; ++j) {
			float repelStep = (float) repel_step;
			if (diffs[j].alpha < exclusion_threshold
				&& diffs[j].L < exclusion_threshold
				&& diffs[j].A < exclusion_threshold
				&& diffs[j].B < exclusion_threshold
				) {

				Lab repel_vec = new Lab();
				repel_vec.alpha = posneg(diffs[j].alpha); repel_vec.L = posneg(diffs[j].L); repel_vec.A = posneg(diffs[j].A); repel_vec.B = posneg(diffs[j].B);

				network[j].alpha += repel_vec.alpha * repelStep;
				network[j].L += repel_vec.L * repelStep;
				network[j].A += repel_vec.A * repelStep;
				network[j].B += repel_vec.B * repelStep;
			}
		}

		repel_points[i] += REPEL_STEP_UP;
	}

	private int contest(int al, float L, float A, float B) {
		/* Calculate the component-wise differences between target_pix colour and every colour in the network, and weight according
		* to component relevance.
		*/
		Lab[] diffs = new Lab[netsize];
		for (int vdx = 0; vdx < netsize; ++vdx) {
			diffs[vdx] = new Lab();
			diffs[vdx].alpha = Math.abs(network[vdx].alpha - al);
			diffs[vdx].L = Math.abs(network[vdx].L - L);
			diffs[vdx].A = Math.abs(network[vdx].A - A);
			diffs[vdx].B = Math.abs(network[vdx].B - B);
		}

		/* finds closest neuron (min dist) and updates freq */
		/* finds best neuron (min dist-bias) and returns position */
		/* for frequently chosen neurons, freq[i] is high and bias[i] is negative */
		/* bias[i] = gamma*((1/netsize)-freq[i]) */

		int bestpos = 0, bestbiaspos = bestpos;
		double bestd = Integer.MAX_VALUE, bestbiasd = bestd;

		/* Using colorimportance(al) here was causing problems with images that were close to monocolor.
		See bug reports: 3149791, 2938728, 2896731 and 2938710
		*/
		boolean perfect = false; /* Is bestpos a perfect match for the colour. */

		for (int i = 0; i < netsize; ++i) {
			if (!perfect) {
				double bestbiasd_biased = bestbiasd + bias[i];
				double dist = 0;

				/*Detect perfect matches. */
				if (diffs[i].alpha < exclusion_threshold
					&& diffs[i].L < exclusion_threshold
					&& diffs[i].A < exclusion_threshold
					&& diffs[i].B < exclusion_threshold
					) {
					perfect = true;

					/* If we don't have a perfect match, the distance between the
						* target and the neuron is the manhattan distance. */
				}
				else {
					dist = diffs[i].alpha;
					dist += diffs[i].L;

					if (dist < bestd || dist < bestbiasd_biased) {
						dist += diffs[i].A;
						dist += diffs[i].B;
					}
				}

				/* See if the current neuron is better. */
				if (dist < bestd) {
					bestd = dist;
					bestpos = i;
				}
				if (dist < bestbiasd_biased) {
					bestbiasd = dist - bias[i];
					bestbiaspos = i;
				}
			}

			/* Age (decay) the current neurons bias and freq values. */
			double betafreq = freq[i] * beta;
			freq[i] -= betafreq;
			bias[i] += betafreq * gamma;
		}

		/* Increase the freq and bias values for the chosen neuron. */
		freq[bestpos] += beta;
		bias[bestpos] -= betagamma;
		
		/* If our bestpos pixel is a 'perfect' match, we return bestpos, not bestbiaspos.  That is, we only decide to look at
		* bestbiaspos if the current target pixel wasn't a good enough match with the bestpos neuron, and there is some hope that
		* we can train the bestbiaspos neuron to become a better match. */
		if (perfect)
			return -bestpos - 1;  /* flag this was a perfect match */

		return bestbiaspos;
	}
	
	protected final void setColorModel(final Color[] palette)
	{
		m_colorModel = BitmapUtilities.createColorModel(palette, m_transparentPixelIndex, hasSemiTransparency, hasSemiTransparency);
	}

	private void learn(final int samplefac, final int nMaxColors, final Color[] pixels) {
		int pos = 0;
		int alphadec = 30 + ((samplefac - 1) / 3);
		final int lengthcount = pixels.length;
		int samplepixels = lengthcount / samplefac;
		int delta = samplepixels / ncycles;  /* here's a problem with small images: samplepixels < ncycles => delta = 0 */
		if (delta == 0)
			delta = 1;        /* kludge to fix */
		double alpha = initalpha;
		double radius = initradius;
		if (netsize <= 2 * initradius)
			radius = netsize / 2;

		int rad = (int) radius;
		if (rad <= 1)
			rad = 0;

		for (int i = 0; i < rad; ++i)
			radpower[i] = Math.floor(alpha * (((BitmapUtilities.sqr(rad) - BitmapUtilities.sqr(i)) * radiusbias) / BitmapUtilities.sqr(rad)));

		int step = (int) (random.nextDouble() * lengthcount);

		int learning_extension = normal_learning_extension_factor;
		if (netsize < extra_long_colour_threshold)
			learning_extension = 2 + ((extra_long_colour_threshold - netsize) / extra_long_divisor);

		int i = 0;
		while (i < learning_extension * samplepixels) {
			Color c = pixels[pos];

			int al = c.getAlpha();
			if (c.getAlpha() <= alphaThreshold)
				c = m_transparentColor;
			Lab lab1 = getLab(c.getRGB());

			int j = contest(al, lab1.L, lab1.A, lab1.B);

			/* Determine if the colour was a perfect match.  j contains a factor encoded boolean. Horrible code to extract it. */
			boolean was_perfect = (j < 0);
			j = (j < 0 ? -(j + 1) : j);

			altersingle(alpha, j, al, lab1.L, lab1.A, lab1.B);
			if (rad > 0 && !was_perfect)
				alterneigh(rad, j, al, lab1.L, lab1.A, lab1.B);   /* alter neighbours */
			else if (rad > 0 && was_perfect)
				repelcoincident(j);  /* repel neighbours in colour space */

			pos += step;
			while (pos >= lengthcount)
				pos -= lengthcount;

			if (++i % delta == 0) {                    /* FPE here if delta=0*/
				alpha -= alpha / (learning_extension * (double) alphadec);
				radius -= radius / (double)radiusdec;
				rad = (int) radius;
				if (rad <= 1)
					rad = 0;
				for (int k = 0; k < rad; ++k)
					radpower[k] = Math.floor(alpha * (((BitmapUtilities.sqr(rad) - BitmapUtilities.sqr(k)) * radiusbias) / BitmapUtilities.sqr(rad)));
			}
		}
	}

	void inxbuild(Color[] palette) {
		int nMaxColors = palette.length;

		int previouscol = 0;
		int startpos = 0;

		for (int i = 0; i < nMaxColors; ++i) {
			int smallpos = i;
			float smallval = network[i].L;			// index on L
											// find smallest in i..netsize-1
			for (int j = i + 1; j < nMaxColors; ++j) {
				if (network[j].L < smallval) {		// index on L				
					smallpos = j;
					smallval = network[j].L;	// index on L
				}
			}
			// swap p (i) and q (smallpos) entries
			if (i != smallpos) {
				Lab temp = network[i];
				network[i] = network[smallpos];
				network[smallpos] = temp;
			}

			// smallval entry is now in position i
			if (smallval != previouscol) {
				netindex[previouscol] = (startpos + i) >> 1;
				for (int j = previouscol + 1; j < smallval; ++j)
					netindex[j] = i;
				previouscol = (int) smallval;
				startpos = i;
			}
		}

		netindex[previouscol] = (startpos + maxnetpos) >> 1;
		for (int j = previouscol + 1; j < netsize; ++j)
			netindex[j] = maxnetpos;

		for (int k = 0; k < nMaxColors; ++k) {
			Lab lab1 = network[k];
			lab1.alpha = round_biased(network[k].alpha);
			palette[k] = CIELABConvertor.LAB2RGB(lab1);
		}
	}

	public short nearestColorIndex(final Color[] palette, Color c, final int pos)
	{
		Short got = nearestMap.get(c.getRGB());
		if (got != null)
			return got;
		
		short k = 0;		
		if (c.getAlpha() <= alphaThreshold)
			c = m_transparentColor;
		if(palette.length > 2 && hasAlpha() && c.getAlpha() > alphaThreshold)
			k = 1;
		
		double mindist = 1e100;
		Lab lab1 = getLab(c.getRGB());
		for (short i=k; i<palette.length; ++i) {
			Color c2 = palette[i];

			double curdist = hasSemiTransparency ? BitmapUtilities.sqr(c2.getAlpha() - c.getAlpha()) / Math.exp(1.5) : 0;
			if (curdist > mindist)
				continue;
			
			Lab lab2 = getLab(c2.getRGB());
			if (palette.length <= 4) {
				curdist = BitmapUtilities.sqr(c2.getRed() - c.getRed())
						+ BitmapUtilities.sqr(c2.getGreen() - c.getGreen()) + BitmapUtilities.sqr(c2.getBlue() - c.getBlue());
				if (hasSemiTransparency)
					curdist += BitmapUtilities.sqr(c2.getAlpha() - c.getAlpha());
			}
			else if (hasSemiTransparency) {		
				curdist += BitmapUtilities.sqr(lab2.L - lab1.L);
				if (curdist > mindist)
					continue;
				
				curdist += BitmapUtilities.sqr(lab2.A - lab1.A);
				if (curdist > mindist)
					continue;
				
				curdist += BitmapUtilities.sqr(lab2.B - lab1.B);
			}
			else if (palette.length > 32) {				
				curdist += Math.abs(lab2.L - lab1.L);
				if (curdist > mindist)
					continue;
				
				curdist += Math.sqrt(BitmapUtilities.sqr(lab2.A - lab1.A) + BitmapUtilities.sqr(lab2.B - lab1.B));
			}
			else {		
				double deltaL_prime_div_k_L_S_L = CIELABConvertor.L_prime_div_k_L_S_L(lab1, lab2);
				curdist += BitmapUtilities.sqr(deltaL_prime_div_k_L_S_L);
				if (curdist > mindist)
					continue;
	
				MutableDouble a1Prime = new MutableDouble(), a2Prime = new MutableDouble(), CPrime1 = new MutableDouble(), CPrime2 = new MutableDouble();
				float deltaC_prime_div_k_L_S_L = CIELABConvertor.C_prime_div_k_L_S_L(lab1, lab2, a1Prime, a2Prime, CPrime1, CPrime2);
				curdist += BitmapUtilities.sqr(deltaC_prime_div_k_L_S_L);
				if (curdist > mindist)
					continue;
	
				MutableDouble barCPrime = new MutableDouble(), barhPrime = new MutableDouble();
				float deltaH_prime_div_k_L_S_L = CIELABConvertor.H_prime_div_k_L_S_L(lab1, lab2, a1Prime, a2Prime, CPrime1, CPrime2, barCPrime, barhPrime);
				curdist +=  BitmapUtilities.sqr(deltaH_prime_div_k_L_S_L);
				if (curdist > mindist)
					continue;
	
				curdist += CIELABConvertor.R_T(barCPrime, barhPrime, deltaC_prime_div_k_L_S_L, deltaH_prime_div_k_L_S_L);
			}

			if (curdist > mindist)
				continue;
			mindist = curdist;
			k = i;
		}
		nearestMap.put(c.getRGB(), k);
		return k;
	}

	void Clear() {
		pixelMap.clear();
		nearestMap.clear();
	}
	
	protected Ditherable getDitherFn(final boolean dither) {
		return new Ditherable() {
			@Override
			public int getColorIndex(Color c) {
				return BitmapUtilities.getColorIndex(c, hasSemiTransparency, hasAlpha());
			}
			
			@Override
			public short nearestColorIndex(Color[] palette, Color c, final int pos) {
				return NeuQuantizer.this.nearestColorIndex(palette, c, pos);
			}
			
		};
	}
	
	protected short[] dither(Color[] palette, int width, int height, boolean dither)
	{
		Ditherable ditherable = getDitherFn(dither);
		if(dither) {
			double weight = 1;
			if(hasSemiTransparency)
				weight *= -1;
			short[] qPixels = GilbertCurve.dither(width, height, pixels, palette, ditherable, null, weight, dither); 
			return qPixels;
		}
		
		short[] qPixels = BitmapUtilities.quantize_image(width, height, pixels, palette, ditherable, hasSemiTransparency, dither);			
		Clear();
		return qPixels;
	}
	
	protected Color[] grabPixels(int[] pixels, int nMaxColors, ThreadLocal<Boolean> isSemiTransparency) {
		final Color[] cPixels = new Color[pixels.length];
		int semiTransCount = 0;
		for (int i = 0; i < pixels.length; ++i) {
			int pixel = pixels[i];
			int alfa = (pixel >> 24) & 0xff;
			cPixels[i] = new Color(pixel, true);
			if (alfa < 0xE0) {	
				if (alfa == 0) {
					m_transparentPixelIndex = i;
					if(nMaxColors > 2)
						m_transparentColor = cPixels[i];
					else
						cPixels[i] = m_transparentColor;
				}
				else if(alfa > alphaThreshold)
					++semiTransCount;
			}		
		}
		
		hasSemiTransparency = semiTransCount > 0;
		if(isSemiTransparency != null)
			isSemiTransparency.set(hasSemiTransparency);
		return cPixels;
	}

	// The work horse for NeuralNet color quantizing.
	public BufferedImage convert(int nMaxColors, boolean dither) {		
		Clear();
		
		Color[] palette = new Color[nMaxColors];
		final Color[] cPixels = grabPixels(pixels, nMaxColors, null);
		if (nMaxColors > 2) {
			netsize = nMaxColors;		// number of colors used
			maxnetpos = netsize - 1;
			initrad = netsize < 8 ? 1 : (netsize >> 3);
			initradius = initrad * 1.0;

			setUpArrays();
			learn(dither ? 5 : 1, nMaxColors, cPixels);
			inxbuild(palette);
		}
		else {
			if (hasAlpha()) {
				palette[0] = m_transparentColor;
				palette[1] = Color.BLACK;
			}
			else {
				palette[0] = Color.BLACK;
				palette[1] = Color.WHITE;
			}
		}

		short[] qPixels = dither(palette, width, height, dither);
		if (hasAlpha() && nMaxColors > 2)
		{
			short k = qPixels[m_transparentPixelIndex];
			palette[k] = m_transparentColor;
		}
		setColorModel(palette);
		return BitmapUtilities.processImagePixels(qPixels, m_colorModel, width, height);
	}

	public boolean hasAlpha() {
		return m_transparentPixelIndex > -1;
	}
}
