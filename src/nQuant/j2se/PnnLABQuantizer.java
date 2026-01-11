package nQuant.j2se;

/* Fast pairwise nearest neighbor based algorithm with CIELAB color space advanced version
Copyright (c) 2018-2026 Miller Cy Chan
* error measure; time used is proportional to number of bins squared - WJ */

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.awt.image.ImageObserver;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;

import nQuant.j2se.CIELABConvertor.Lab;
import nQuant.j2se.CIELABConvertor.MutableDouble;

public class PnnLABQuantizer extends PnnQuantizer {
	protected float[] saliencies;
	protected Map<Integer, Lab> pixelMap = new HashMap<>();
	
	private static Random random = new Random();
	private Color[] cPixels;
	private boolean isGA = false;
	private double proportional, ratioY = .5;

	public PnnLABQuantizer(BufferedImage im, ImageObserver obs) {
		super(im, obs);
	}
	
	PnnLABQuantizer(PnnLABQuantizer quantizer) {
		super(quantizer);
		saliencies = quantizer.saliencies;
		pixelMap = new HashMap<>(quantizer.pixelMap);
		cPixels = quantizer.cPixels;
		isGA = true;
		proportional = quantizer.proportional;
	}

	private static final class Pnnbin {
		float ac = 0, Lc = 0, Ac = 0, Bc = 0;
		float cnt = 0, err = 0;
		int nn, fw, bk, tm, mtm;
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

	private void find_nn(Pnnbin[] bins, int idx, boolean texicab)
	{
		int nn = 0;
		double err = 1e100;

		Pnnbin bin1 = bins[idx];
		float n1 = bin1.cnt;
		
		Lab lab1 = new Lab();
		lab1.alpha = bin1.ac; lab1.L = bin1.Lc; lab1.A = bin1.Ac; lab1.B = bin1.Bc;
		for (int i = bin1.fw; i != 0; i = bins[i].fw) {
			float n2 = bins[i].cnt;
			double nerr2 = (n1 * n2) / (n1 + n2);
			if (nerr2 >= err)
				continue;
			
			Lab lab2 = new Lab();
			lab2.alpha = bins[i].ac; lab2.L = bins[i].Lc; lab2.A = bins[i].Ac; lab2.B = bins[i].Bc;
			double alphaDiff = hasSemiTransparency ? BitmapUtilities.sqr(lab2.alpha - lab1.alpha) / Math.exp(1.75) : 0;
			double nerr = nerr2 * alphaDiff;
			if (nerr >= err)
				continue;
			
			if(!texicab) {
				nerr += (1 - ratio) * nerr2 * BitmapUtilities.sqr(lab2.L - lab1.L);
				if (nerr >= err)
					continue;

				nerr += (1 - ratio) * nerr2 * BitmapUtilities.sqr(lab2.A - lab1.A);
				if (nerr >= err)
					continue;

				nerr += (1 - ratio) * nerr2 * BitmapUtilities.sqr(lab2.B - lab1.B);
			}
			else {
				nerr += (1 - ratio) * nerr2 * Math.abs(lab2.L - lab1.L);
				if (nerr >= err)
					continue;

				nerr += (1 - ratio) * nerr2 * Math.sqrt(BitmapUtilities.sqr(lab2.A - lab1.A) + BitmapUtilities.sqr(lab2.B - lab1.B));
			}

			if (nerr > err)
				continue;
			
			float deltaL_prime_div_k_L_S_L = CIELABConvertor.L_prime_div_k_L_S_L(lab1, lab2);
			nerr += ratio * nerr2 * BitmapUtilities.sqr(deltaL_prime_div_k_L_S_L);
			if (nerr > err)
				continue;

			MutableDouble a1Prime = new MutableDouble(), a2Prime = new MutableDouble(), CPrime1 = new MutableDouble(), CPrime2 = new MutableDouble();
			float deltaC_prime_div_k_L_S_L = CIELABConvertor.C_prime_div_k_L_S_L(lab1, lab2, a1Prime, a2Prime, CPrime1, CPrime2);
			nerr += ratio * nerr2 * BitmapUtilities.sqr(deltaC_prime_div_k_L_S_L);
			if (nerr > err)
				continue;

			MutableDouble barCPrime = new MutableDouble(), barhPrime = new MutableDouble();
			float deltaH_prime_div_k_L_S_L = CIELABConvertor.H_prime_div_k_L_S_L(lab1, lab2, a1Prime, a2Prime, CPrime1, CPrime2, barCPrime, barhPrime);
			nerr += ratio * nerr2 * BitmapUtilities.sqr(deltaH_prime_div_k_L_S_L);
			if (nerr > err)
				continue;

			nerr += ratio * nerr2 * CIELABConvertor.R_T(barCPrime, barhPrime, deltaC_prime_div_k_L_S_L, deltaH_prime_div_k_L_S_L);
			if (nerr > err)
				continue;

			err = nerr;
			nn = i;
		}
		bin1.err = (float) err;
		bin1.nn = nn;
	}

	protected QuanFn getQuanFn(int nMaxColors, short quan_rt) {
		if (quan_rt > 0) {
			if (quan_rt > 1)
				return cnt -> (float) Math.pow(cnt, 0.75);
			if (nMaxColors < 64)
				return cnt -> (int) Math.sqrt(cnt);

			return cnt -> (float) Math.sqrt(cnt);
		}
		return cnt -> cnt;
	}

	@Override
	protected Color[] pnnquan(final Color[] pixels, int nMaxColors)
	{
		short quan_rt = (short) 1;
		Pnnbin[] bins = new Pnnbin[65536];
		saliencies = nMaxColors >= 128 ? null : new float[pixels.length];
		float saliencyBase = .1f;

		/* Build histogram */
		for (int i = 0; i < pixels.length; ++i) {
			Color c = pixels[i];
			if(c.getAlpha() <= alphaThreshold)
				c = m_transparentColor;
			
			int index = BitmapUtilities.getColorIndex(c, hasSemiTransparency, hasAlpha());
			Lab lab1 = getLab(c.getRGB());
			
			if(bins[index] == null)
				bins[index] = new Pnnbin();
			Pnnbin tb = bins[index];
			tb.ac += lab1.alpha;
			tb.Lc += lab1.L;
			tb.Ac += lab1.A;
			tb.Bc += lab1.B;
			tb.cnt += 1.0f;
			if(saliencies != null)
				saliencies[i] = saliencyBase + (1 - saliencyBase) * lab1.L / 100f * lab1.alpha / 255f;
		}

		/* Cluster nonempty bins at one end of array */
		int maxbins = 0;

		for (int i = 0; i < bins.length; ++i) {
			if (bins[i] == null)
				continue;

			float d = 1f / bins[i].cnt;
			bins[i].ac *= d;
			bins[i].Lc *= d;
			bins[i].Ac *= d;
			bins[i].Bc *= d;

			bins[maxbins++] = bins[i];
		}
		
		proportional = BitmapUtilities.sqr(nMaxColors) / maxbins;
		if ((hasAlpha() || hasSemiTransparency) && nMaxColors < 32)
			quan_rt = -1;
		
		weight = Math.min(0.9, nMaxColors * 1.0 / maxbins);
		if ((nMaxColors < 16 && weight < .0075) || weight < .001 || (weight > .0015 && weight < .0022))
			quan_rt = 2;
		if (weight < (isGA ? .03 : .04) && PG < 1 && PG >= coeffs[0][1]) {
			if (nMaxColors >= 64)
				quan_rt = 0;
		}
		if (nMaxColors > 16 && nMaxColors < 64) {
			double weightB = nMaxColors / 8000.0;
			if (Math.abs(weightB - weight) < .001)
				quan_rt = 2;
		}

		if(pixelMap.size() <= nMaxColors) {
			/* Fill palette */
			Color[] palette = new Color[pixelMap.size()];
			int k = 0;
			for (Integer pixel : pixelMap.keySet()) {
				Color c = new Color(pixel, hasAlpha());
				palette[k++] = c;

				if(k > 1 && c.getAlpha() == 0) {
					palette[k - 1] = palette[0]; palette[0] = c;
				}
			}

			return palette;
		}
		
		QuanFn quanFn = getQuanFn(nMaxColors, quan_rt);
		
		int j = 0;
		for (; j < maxbins - 1; ++j) {
			bins[j].fw = j + 1;
			bins[j + 1].bk = j;
			
			bins[j].cnt = quanFn.get(bins[j].cnt);
		}
		bins[j].cnt = quanFn.get(bins[j].cnt);

		final boolean texicab = proportional > .0225 && !hasSemiTransparency;

		if(!isGA) {
			if(hasSemiTransparency)
				ratio = .5;
			else if(quan_rt != 0 && nMaxColors < 64) {
				if (proportional > .018 && proportional < .022)
					ratio = Math.min(1.0, proportional + weight * Math.exp(3.13));
				else if(proportional > .1)
					ratio = Math.min(1.0, 1.0 - weight);
				else if(proportional > .04)
					ratio = Math.min(1.0, weight * Math.exp(1.56));
				else if(proportional > .025 && (weight < .002 || weight > .0022))
					ratio = Math.min(1.0, proportional + weight * Math.exp(3.66));
				else
					ratio = Math.min(1.0, proportional + weight * Math.exp(1.718));
			}
			else if(nMaxColors > 256)
				ratio = Math.min(1.0, 1 - 1.0 / proportional);
			else
				ratio = Math.min(1.0, 1 - weight * .7);
			
			if (!hasSemiTransparency && quan_rt < 0)
				ratio = Math.min(1.0, weight * Math.exp(3.13));
		}

		int h, l, l2;
		/* Initialize nearest neighbors and build heap of them */
		int[] heap = new int[bins.length + 1];
		for (int i = 0; i < maxbins; ++i) {
			find_nn(bins, i, texicab);
			/* Push slot on heap */
			float err = bins[i].err;
			for (l = ++heap[0]; l > 1; l = l2) {
				l2 = l >> 1;
				if (bins[h = heap[l2]].err <= err)
					break;
				heap[l] = h;
			}
			heap[l] = i;
		}

		if (!isGA && quan_rt > 0 && nMaxColors < 64 && proportional > .035 && proportional < .1) {
			final int dir = proportional > .04 ? 1 : -1;
			final double margin = dir > 0 ? .002 : .0025;
			final double delta = weight > margin && weight < .003 ? 1.872 : 1.632;
			ratio = Math.min(1.0, proportional + dir * weight * Math.exp(delta));
		}
		else if (isGA)
			ratio = ratioY;

		/* Merge bins which increase error the least */
		int extbins = maxbins - nMaxColors;
		for (int i = 0; i < extbins; ) {
			Pnnbin tb = null;
			/* Use heap to find which bins to merge */
			for (;;) {
				int b1 = heap[1];
				tb = bins[b1]; /* One with least error */
				/* Is stored error up to date? */
				if ((tb.tm >= tb.mtm) && (bins[tb.nn].mtm <= tb.tm))
					break;
				if (tb.mtm == 0xFFFF) /* Deleted node */
					b1 = heap[1] = heap[heap[0]--];
				else /* Too old error value */
				{
					find_nn(bins, b1, texicab);
					tb.tm = i;
				}
				/* Push slot down */
				float err = bins[b1].err;
				for (l = 1; (l2 = l + l) <= heap[0]; l = l2) {
					if ((l2 < heap[0]) && (bins[heap[l2]].err > bins[heap[l2 + 1]].err))
						++l2;
					if (err <= bins[h = heap[l2]].err)
						break;
					heap[l] = h;
				}
				heap[l] = b1;
			}

			/* Do a merge */
			Pnnbin nb = bins[tb.nn];
			float n1 = tb.cnt;
			float n2 = nb.cnt;
			float d = 1.0f / (n1 + n2);
			tb.ac = d * (n1 * tb.ac + n2 * nb.ac);
			tb.Lc = d * (n1 * tb.Lc + n2 * nb.Lc);
			tb.Ac = d * (n1 * tb.Ac + n2 * nb.Ac);
			tb.Bc = d * (n1 * tb.Bc + n2 * nb.Bc);
			tb.cnt += n2;
			tb.mtm = ++i;

			/* Unchain deleted bin */
			bins[nb.bk].fw = nb.fw;
			bins[nb.fw].bk = nb.bk;
			nb.mtm = 0xFFFF;
		}

		/* Fill palette */
		Color[] palette = new Color[extbins > 0 ? nMaxColors : maxbins];
		int k = 0;
		for (int i = 0; k < palette.length; ++k) {
			Lab lab1 = new Lab();
			lab1.alpha = (int) Math.rint(bins[i].ac);
			lab1.L = bins[i].Lc; lab1.A = bins[i].Ac; lab1.B = bins[i].Bc;
			palette[k] = CIELABConvertor.LAB2RGB(lab1);

			i = bins[i].fw;
		}

		return palette;
	}

	@Override
	public short nearestColorIndex(final Color[] palette, Color c, final int pos)
	{
		final int offset = getColorIndex(c);
		int got = nearestMap[offset];
		if (got > 0)
			return (short) (got - 1);

		short k = 0;
		if (c.getAlpha() <= alphaThreshold)
			c = m_transparentColor;
		if(palette.length > 2 && hasAlpha() && c.getAlpha() > alphaThreshold)
			k = 1;
		
		double mindist = 1e100;
		Lab lab1 = getLab(c.getRGB());
		for (short i = k; i < palette.length; ++i) {
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
			else if (hasSemiTransparency || palette.length < 16) {
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
		nearestMap[offset] = k + 1;
		return k;
	}
	
	protected short hybridColorIndex(final Color[] palette, Color c, final int pos)
	{
		final int offset = getColorIndex(c);
		int got = nearestMap[offset];
		if (got > 0)
			return (short) (got - 1);

		short k = 0;
		
		double mindist = 1e100;
		Lab lab1 = getLab(c.getRGB());
		for (short i = k; i < palette.length; ++i) {
			Color c2 = palette[i];

			double curdist = hasSemiTransparency ? BitmapUtilities.sqr(c2.getAlpha() - c.getAlpha()) / Math.exp(1.5) : 0;
			
			Lab lab2 = getLab(c2.getRGB());
			if (Math.abs(lab2.L - lab1.L) < palette.length || saliencies[pos] < .2 || saliencies[pos] > .8) {
				curdist += BitmapUtilities.sqr(lab2.L - lab1.L);
				if (curdist > mindist)
					continue;
				
				curdist += BitmapUtilities.sqr(lab2.A - lab1.A);
				if (curdist > mindist)
					continue;
				
				curdist += BitmapUtilities.sqr(lab2.B - lab1.B);
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
		nearestMap[offset] = k + 1;
		return k;
	}

	@Override
	protected short closestColorIndex(final Color[] palette, Color c, final int pos)
	{
		if(PG < coeffs[0][1] && BlueNoise.TELL_BLUE_NOISE[pos & 4095] > -88)
			return hybridColorIndex(palette, c, pos);
		
		if (c.getAlpha() <= alphaThreshold)
			return nearestColorIndex(palette, c, pos);
		
		final int offset = getColorIndex(c);
		int[] closest = closestMap.get(offset);
		if (closest == null) {
			closest = new int[4];
			closest[2] = closest[3] = Integer.MAX_VALUE;

			for (short k = 0; k < palette.length; ++k) {
				Color c2 = palette[k];
				
				double err = PR * (1 - ratio) * BitmapUtilities.sqr(c2.getRed() - c.getRed());
				if (err >= closest[3])
					continue;

				err += PG * (1 - ratio) * BitmapUtilities.sqr(c2.getGreen() - c.getGreen());
				if (err >= closest[3])
					continue;
				
				err += PB * (1 - ratio) * BitmapUtilities.sqr(c2.getBlue() - c.getBlue());
				if (err >= closest[3])
					continue;

				if(hasSemiTransparency)
					err += PA * BitmapUtilities.sqr(c2.getAlpha() - c.getAlpha());

				for (int i = 0; i < coeffs.length; ++i) {
					err += ratio * BitmapUtilities.sqr(coeffs[i][0] * (c2.getRed() - c.getRed()));
					if (err >= closest[3])
						break;
					
					err += ratio * BitmapUtilities.sqr(coeffs[i][1] * (c2.getGreen() - c.getGreen()));
					if (err >= closest[3])
						break;
					
					err += ratio * BitmapUtilities.sqr(coeffs[i][2] * (c2.getBlue() - c.getBlue()));
					if (err >= closest[3])
						break;
				}

				if (err < closest[2]) {
					closest[1] = closest[0];
					closest[3] = closest[2];
					closest[0] = k;
					closest[2] = (int) err;
				}
				else if (err < closest[3]) {
					closest[1] = k;
					closest[3] = (int) err;
				}
			}

			if (closest[3] == Integer.MAX_VALUE)
				closest[1] = closest[0];

			closestMap.put(offset, closest);
		}
		
		int idx = 1;
		if (closest[2] == 0 || (random.nextInt(closest[3] + closest[2])) <= closest[3])
			idx = 0;
		
		int MAX_ERR = palette.length;
		if(closest[idx + 2] >= MAX_ERR || closest[idx] == 0 || palette[closest[idx]].getAlpha() < c.getAlpha())
			return nearestColorIndex(palette, c, pos);
		return (short) closest[idx];
	}	

	protected Ditherable getDitherFn() {
		return new Ditherable() {
			@Override
			public int getColorIndex(Color c) {
				return BitmapUtilities.getColorIndex(c, hasSemiTransparency, hasAlpha());
			}

			@Override
			public short nearestColorIndex(Color[] palette, Color c, final int pos) {
				if (palette.length <= 4)
					return PnnLABQuantizer.this.nearestColorIndex(palette, c, pos);
				if (isGA() && palette.length < 16)
					return PnnLABQuantizer.this.hybridColorIndex(palette, c, pos);
				return PnnLABQuantizer.this.nearestColorIndex(palette, c, pos);
			}
			
		};
	}

	protected void clear()
	{
		m_palette = null;
		saliencies = null;
		closestMap.clear();
		nearestMap = new int[65536];
	}

	@Override
	protected short[] dither(Color[] palette, int width, int height, boolean dither)
	{
		Ditherable ditherable = getDitherFn();
		if(hasSemiTransparency)
			weight *= -1;
		
		if(dither && saliencies == null && (palette.length <= 256 || weight > .99)) {
			saliencies = new float[pixels.length];
			float saliencyBase = .1f;

			for (int i = 0; i < pixels.length; ++i) {
				Color c = cPixels[i];
				Lab lab1 = getLab(c.getRGB());

				saliencies[i] = saliencyBase + (1 - saliencyBase) * lab1.L / 100f * lab1.alpha / 255f;
			}
		}
		short[] qPixels = GilbertCurve.dither(width, height, pixels, palette, ditherable, saliencies, weight, dither);

		if(!dither && palette.length > 32) {
			double delta = BitmapUtilities.sqr(palette.length) / pixelMap.size();
			float weight = delta > 0.023 ? 1.0f : (float) (37.013 * delta + 0.906);
			BlueNoise.dither(width, height, pixels, palette, ditherable, qPixels, weight);
		}
		
		return qPixels;
	}

	@Override
	protected Color[] grabPixels(int[] pixels, int nMaxColors, ThreadLocal<Boolean> isSemiTransparency) {
		cPixels = super.grabPixels(pixels, nMaxColors, isSemiTransparency);
		return cPixels;
	}

	int[] getPixels()
	{
		return pixels;
	}

	void setPixels(Color[] cPixels)
	{
		pixels = new int[cPixels.length];
		int i = 0;
		for(Color c : cPixels)
			pixels[i++] = c.getRGB();
	}

	Color[] getCPixels()
	{
		return cPixels;
	}
	
	boolean isGA() {
		return isGA;
	}

	void setRatio(double ratioX, double ratioY) {
		this.ratio = Math.min(1.0, ratioX);
		this.ratioY = Math.min(1.0, ratioY);
		clear();
	}
	
}
