package nQuant.j2se;

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.awt.image.ImageObserver;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;

import nQuant.j2se.CIELABConvertor.Lab;
import nQuant.j2se.CIELABConvertor.MutableDouble;

public class PnnLABQuantizer extends PnnQuantizer {
	private double PR = .2126, PG = .7152, PB = .0722;
	private double ratio = 1.0;
	private Map<Integer, Lab> pixelMap = new HashMap<Integer, Lab>();	

	public PnnLABQuantizer(BufferedImage im, ImageObserver obs) throws IOException {
		super(im, obs);
	}

	private static final class Pnnbin {
		float ac = 0, Lc = 0, Ac = 0, Bc = 0;
		float cnt = 0, err = 0;
		int nn, fw, bk, tm, mtm;
	}

	private Lab getLab(final int pixel)
	{
		Lab lab1 = pixelMap.get(pixel);
		if (lab1 == null) {
			lab1 = CIELABConvertor.RGB2LAB(pixel);
			pixelMap.put(pixel, lab1);
		}
		return lab1;
	}

	private void find_nn(Pnnbin[] bins, int idx, int nMaxColors)
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
			double alphaDiff = hasSemiTransparency ? Math.abs(lab2.alpha - lab1.alpha) : 0;
			double nerr = nerr2 * BitmapUtilities.sqr(alphaDiff) / Math.exp(1.5);
			if (nerr >= err)
				continue;
			
			nerr += (1 - ratio) * nerr2 * BitmapUtilities.sqr(lab2.L - lab1.L);
			if (nerr >= err)
				continue;

			nerr += (1 - ratio) * nerr2 * BitmapUtilities.sqr(lab2.A - lab1.A);
			if (nerr >= err)
				continue;

			nerr += (1 - ratio) * nerr2 * BitmapUtilities.sqr(lab2.B - lab1.B);

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

	@Override
	protected Color[] pnnquan(final Color[] pixels, int nMaxColors, short quan_rt)
	{
		if(hasSemiTransparency)
			PR = PG = PB = 1.0;
		else if(width < 512 || height < 512) {
			PR = 0.299; PG = 0.587; PB = 0.114;
		}
		
		Pnnbin[] bins = new Pnnbin[65536];		

		/* Build histogram */
		for (Color c : pixels) {
			// !!! Can throw gamma correction in here, but what to do about perceptual
			// !!! nonuniformity then?
			
			int index = BitmapUtilities.getColorIndex(c, hasSemiTransparency, m_transparentPixelIndex >= 0);
			Lab lab1 = getLab(c.getRGB());
			if(bins[index] == null)
				bins[index] = new Pnnbin();
			Pnnbin tb = bins[index];
			tb.ac += c.getAlpha();
			tb.Lc += lab1.L;
			tb.Ac += lab1.A;
			tb.Bc += lab1.B;
			tb.cnt += 1.0f;
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
		
		double proportional = BitmapUtilities.sqr(nMaxColors) / maxbins;
		if ((m_transparentPixelIndex >= 0 || hasSemiTransparency) && nMaxColors < 32)
			quan_rt = -1;
		else if ((proportional < .018 || proportional > .5) && nMaxColors < 64)
			quan_rt = 0;
		
		int j = 0;
		for (; j < maxbins - 1; ++j) {
			bins[j].fw = j + 1;
			bins[j + 1].bk = j;
			
			if (quan_rt > 0) {
				bins[j].cnt = (float) Math.sqrt(bins[j].cnt);
				if (nMaxColors < 64)
					bins[j].cnt = (int)bins[j].cnt;
			}
		}
		if (quan_rt > 0) {
			bins[j].cnt = (float) Math.sqrt(bins[j].cnt);
			if (nMaxColors < 64)
				bins[j].cnt = (int)bins[j].cnt;
		}
		
		if(quan_rt != 0 && nMaxColors < 64) {
			if (proportional > .018 && proportional < .022)
				ratio = Math.min(1.0, proportional + nMaxColors * Math.exp(3.13) / maxbins);
			else if(proportional > .1)
				ratio = Math.min(1.0, proportional + nMaxColors * Math.exp(2.68) / maxbins);
			else
				ratio = Math.min(1.0, proportional + nMaxColors * Math.exp(1.718) / maxbins);
		}
		else
			ratio = Math.min(1.0, 0.14 * Math.exp(4.679 * proportional));
		
		if (quan_rt < 0) {
			ratio += 0.5;		
			ratio = Math.min(1.0, ratio);
		}
		
		int h, l, l2;
		/* Initialize nearest neighbors and build heap of them */
		int[] heap = new int[bins.length + 1];
		for (int i = 0; i < maxbins; ++i) {
			find_nn(bins, i, nMaxColors);
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

		if (quan_rt > 0 && nMaxColors < 64 && proportional > .035) {
			int dir = proportional > .04 ? 1 : -1;
			ratio = Math.min(1.0, proportional + dir * nMaxColors * Math.exp(1.632) / maxbins);
		}

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
					find_nn(bins, b1, nMaxColors);
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
			tb.cnt += nb.cnt;
			tb.mtm = ++i;

			/* Unchain deleted bin */
			bins[nb.bk].fw = nb.fw;
			bins[nb.fw].bk = nb.bk;
			nb.mtm = 0xFFFF;
		}

		/* Fill palette */
		Color[] palette = new Color[nMaxColors];
		int k = 0;
		for (int i = 0; ; ++k) {
			Lab lab1 = new Lab();
			lab1.alpha = (int) Math.rint(bins[i].ac);
			lab1.L = bins[i].Lc; lab1.A = bins[i].Ac; lab1.B = bins[i].Bc;
			palette[k] = CIELABConvertor.LAB2RGB(lab1);
			if (m_transparentPixelIndex >= 0 && lab1.alpha == 0) {
				Color temp = palette[0]; palette[0] = m_transparentColor; palette[k] = temp;
			}

			if ((i = bins[i].fw) == 0)
				break;
		}

		setColorModel(palette);		
		return palette;
	}

	@Override
	public short nearestColorIndex(final Color[] palette, final Color c)
	{
		Short got = nearestMap.get(c.getRGB());
		if (got != null)
			return got;
		
		short k = 0;
		if (c.getAlpha() <= alphaThreshold)
            return k;
		
		double mindist = 1e100;
		Lab lab1 = getLab(c.getRGB());
		for (short i=0; i<palette.length; ++i) {
			Color c2 = palette[i];
			if(c2 == null)
				break;

			double curdist = hasSemiTransparency ? BitmapUtilities.sqr(c2.getAlpha() - c.getAlpha()) : 0;
			if (curdist > mindist)
				continue;
			
			if (palette.length > 32 || palette.length <= 4 || hasSemiTransparency) {
				curdist += PR * BitmapUtilities.sqr(c2.getRed() - c.getRed());
				if (curdist > mindist)
					continue;

				curdist += PG * BitmapUtilities.sqr(c2.getGreen() - c.getGreen());
				if (curdist > mindist)
					continue;

				curdist += PB * BitmapUtilities.sqr(c2.getBlue() - c.getBlue());
				if (PB < 1) {
					if (curdist > mindist)
						continue;

					Lab lab2 = getLab(c2.getRGB());
					curdist += BitmapUtilities.sqr(lab2.B - lab1.B) / 2.0;
				}
			}
			else {				
				Lab lab2 = getLab(c2.getRGB());
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

	@Override
	protected short closestColorIndex(final Color[] palette, final Color c)
	{
		int[] closest = closestMap.get(c.getRGB());
		if (closest == null) {
			closest = new int[4];
			closest[2] = closest[3] = Short.MAX_VALUE;

			for (short k = 0; k < palette.length; ++k) {
				Color c2 = palette[k];
				if(c2 == null)
					break;
				
				double err = PR * BitmapUtilities.sqr(c2.getRed() - c.getRed()) + PG * BitmapUtilities.sqr(c2.getGreen() - c.getGreen()) + PB * BitmapUtilities.sqr(c2.getBlue() - c.getBlue());				
				if (err < closest[2]) {
					closest[1] = closest[0];
					closest[3] = closest[2];
					closest[0] = k;
					closest[2] = (short) err;
				}
				else if (err < closest[3]) {
					closest[1] = k;
					closest[3] = (short) err;
				}
			}

			if (closest[3] == Short.MAX_VALUE)
				closest[1] = closest[0];
			
			closestMap.put(c.getRGB(), closest);
		}

		int MAX_ERR = palette.length;
		Random rand = new Random();
		if (closest[2] == 0 || (rand.nextInt(32767) % (closest[3] + closest[2])) <= closest[3]) {
			if(closest[2] > MAX_ERR)
				return nearestColorIndex(palette, c);
			return (short) closest[0];
		}
		
		if(closest[3] > MAX_ERR)
			return nearestColorIndex(palette, c);
		return (short) closest[1];
	}	
	
	protected Ditherable getDitherFn() {
		return new Ditherable() {
			@Override
			public int getColorIndex(Color c) {
				return BitmapUtilities.getColorIndex(c, hasSemiTransparency, m_transparentPixelIndex >= 0);
			}
			
			@Override
			public short nearestColorIndex(Color[] palette, Color c) {
				if(hasSemiTransparency)
					return PnnLABQuantizer.this.nearestColorIndex(palette, c);
				return PnnLABQuantizer.this.closestColorIndex(palette, c);
			}
			
		};
	}
	
	@Override
	protected short[] dither(final Color[] cPixels, Color[] palette, int nMaxColors, int width, int height, boolean dither)
    {		
		short[] qPixels;
		if(hasSemiTransparency)
			qPixels = GilbertCurve.dither(width, height, cPixels, palette, getDitherFn(dither), 1.75f);
		else if (nMaxColors < 64 && nMaxColors > 32)
			qPixels = BitmapUtilities.quantize_image(width, height, cPixels, palette, getDitherFn(dither), hasSemiTransparency, dither);
		else if(nMaxColors <= 32)
			qPixels = GilbertCurve.dither(width, height, cPixels, palette, getDitherFn(), 1.5f); 
		else			
			qPixels = GilbertCurve.dither(width, height, cPixels, palette, getDitherFn());			
		
		if(!dither) {
			double delta = BitmapUtilities.sqr(nMaxColors) / pixelMap.size();
			float weight = delta > 0.023 ? 1.0f : (float) (37.013 * delta + 0.906);
			BlueNoise.dither(width, height, cPixels, palette, getDitherFn(), qPixels, weight);
		}
		
		closestMap.clear();
		nearestMap.clear();
		pixelMap.clear();
		return qPixels;
    }

}