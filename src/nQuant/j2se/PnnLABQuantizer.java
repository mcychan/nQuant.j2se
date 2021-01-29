package nQuant.j2se;

import java.awt.Color;
import java.awt.Image;
import java.awt.image.ImageObserver;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import nQuant.j2se.CIELABConvertor.Lab;
import nQuant.j2se.CIELABConvertor.MutableDouble;

public class PnnLABQuantizer extends PnnQuantizer {
	private double PR = .2126, PG = .7152, PB = .0722;
	private double ratio = 1.0;
	private Map<Color, Lab> pixelMap = new HashMap<Color, Lab>();	

	public PnnLABQuantizer(Image im, int w, int h) throws IOException {
		super(im, w, h);
	}

	public PnnLABQuantizer(Image im, ImageObserver obs) throws IOException {
		super(im, obs);
	}

	private static final class Pnnbin {
		double ac = 0, Lc = 0, Ac = 0, Bc = 0, err = 0;
		int cnt = 0;
		int nn, fw, bk, tm, mtm;
	}

	private Lab getLab(final Color c)
	{
		Lab lab1 = pixelMap.get(c);
		if (lab1 == null) {
			lab1 = CIELABConvertor.RGB2LAB(c);
			pixelMap.put(c, lab1);
		}
		return lab1;
	}

	private void find_nn(Pnnbin[] bins, int idx, int nMaxColors)
	{
		int nn = 0;
		double err = Integer.MAX_VALUE;

		Pnnbin bin1 = bins[idx];
		int n1 = bin1.cnt;
		
		Lab lab1 = new Lab();
		lab1.alpha = bin1.ac; lab1.L = bin1.Lc; lab1.A = bin1.Ac; lab1.B = bin1.Bc;
		boolean crossover = Math.random() < ratio;
		for (int i = bin1.fw; i != 0; i = bins[i].fw) {
			double n2 = bins[i].cnt;
			double nerr2 = (n1 * n2) / (n1 + n2);
			if (nerr2 >= err)
				continue;
			
			Lab lab2 = new Lab();
			lab2.alpha = bins[i].ac; lab2.L = bins[i].Lc; lab2.A = bins[i].Ac; lab2.B = bins[i].Bc;
			double alphaDiff = lab2.alpha - lab1.alpha;
			double nerr = nerr2 * sqr(alphaDiff) * alphaDiff / 3.0;
			if (nerr >= err)
				continue;
			
			if (crossover) {
				double deltaL_prime_div_k_L_S_L = CIELABConvertor.L_prime_div_k_L_S_L(lab1, lab2);
				nerr += nerr2 * sqr(deltaL_prime_div_k_L_S_L);
				if (nerr > err)
					continue;
	
				MutableDouble a1Prime = new MutableDouble(), a2Prime = new MutableDouble(), CPrime1 = new MutableDouble(), CPrime2 = new MutableDouble();
				double deltaC_prime_div_k_L_S_L = CIELABConvertor.C_prime_div_k_L_S_L(lab1, lab2, a1Prime, a2Prime, CPrime1, CPrime2);
				nerr += nerr2 * sqr(deltaC_prime_div_k_L_S_L);
				if (nerr > err)
					continue;
	
				MutableDouble barCPrime = new MutableDouble(), barhPrime = new MutableDouble();
				double deltaH_prime_div_k_L_S_L = CIELABConvertor.H_prime_div_k_L_S_L(lab1, lab2, a1Prime, a2Prime, CPrime1, CPrime2, barCPrime, barhPrime);
				nerr += nerr2 * sqr(deltaH_prime_div_k_L_S_L);
				if (nerr > err)
					continue;
	
				nerr += nerr2 * CIELABConvertor.R_T(barCPrime, barhPrime, deltaC_prime_div_k_L_S_L, deltaH_prime_div_k_L_S_L);
			}
			else {
				nerr += nerr2 * sqr(lab2.L - lab1.L);
				if (nerr >= err)
					continue;

				nerr += nerr2 * sqr(lab2.A - lab1.A);
				if (nerr >= err)
					continue;

				nerr += nerr2 * sqr(lab2.B - lab1.B);
			}
			
			if (nerr > err)
				continue;

			err = nerr;
			nn = i;
		}
		bin1.err = err;
		bin1.nn = nn;
	}

	private Color[] pnnquan(final Color[] pixels, int nMaxColors, boolean quan_sqrt)
	{
		Pnnbin[] bins = new Pnnbin[65536];
		int[] heap = new int[65537];
		double err, n1, n2;

		/* Build histogram */
		for (final Color pixel : pixels) {
			// !!! Can throw gamma correction in here, but what to do about perceptual
			// !!! nonuniformity then?
			int index = getColorIndex(pixel, hasSemiTransparency);
			Lab lab1 = getLab(pixel);
			if(bins[index] == null)
				bins[index] = new Pnnbin();
			Pnnbin tb = bins[index];
			tb.ac += pixel.getAlpha();
			tb.Lc += lab1.L;
			tb.Ac += lab1.A;
			tb.Bc += lab1.B;
			tb.cnt++;
		}

		/* Cluster nonempty bins at one end of array */
		int maxbins = 0;

		for (int i = 0; i < bins.length; ++i) {
			if (bins[i] == null)
				continue;

			double d = 1.0 / (double)bins[i].cnt;
			bins[i].ac *= d;
			bins[i].Lc *= d;
			bins[i].Ac *= d;
			bins[i].Bc *= d;
			
			if(quan_sqrt)
				bins[i].cnt = (int) Math.sqrt(bins[i].cnt);
			bins[maxbins++] = bins[i];
		}

		for (int i = 0; i < maxbins - 1; i++) {
			bins[i].fw = (i + 1);
			bins[i + 1].bk = i;
		}
		// !!! Already zeroed out by calloc()
		//	bins[0].bk = bins[i].fw = 0;

		int h, l, l2;
		ratio = sqr(nMaxColors) / pixelMap.size();
		/* Initialize nearest neighbors and build heap of them */
		for (int i = 0; i < maxbins; i++) {
			find_nn(bins, i, nMaxColors);
			/* Push slot on heap */
			err = bins[i].err;
			for (l = ++heap[0]; l > 1; l = l2) {
				l2 = l >> 1;
				if (bins[h = heap[l2]].err <= err)
					break;
				heap[l] = h;
			}
			heap[l] = i;
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
				err = bins[b1].err;
				for (l = 1; (l2 = l + l) <= heap[0]; l = l2) {
					if ((l2 < heap[0]) && (bins[heap[l2]].err > bins[heap[l2 + 1]].err))
						l2++;
					if (err <= bins[h = heap[l2]].err)
						break;
					heap[l] = h;
				}
				heap[l] = b1;
			}

			/* Do a merge */
			Pnnbin nb = bins[tb.nn];
			n1 = tb.cnt;
			n2 = nb.cnt;
			double d = 1.0 / (n1 + n2);
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
		List<Color> palette = new ArrayList<Color>();
		int k = 0;
		for (int i = 0;; ++k) {
			Lab lab1 = new Lab();
			lab1.alpha = (int) Math.rint(bins[i].ac);
			lab1.L = bins[i].Lc; lab1.A = bins[i].Ac; lab1.B = bins[i].Bc;
			palette.add(CIELABConvertor.LAB2RGB(lab1));
			if (m_transparentPixelIndex >= 0 && palette.get(k).equals(m_transparentColor))
				Collections.swap(palette, 0, k);

			if ((i = bins[i].fw) == 0)
				break;
		}

		setColorModel(palette);
		m_palette = palette.toArray(new Color[0]);
		return m_palette;
	}

	private short nearestColorIndex(final Color[] palette, final Color c)
	{
		short k = 0;
		double mindist = SHORT_MAX;
		Lab lab1 = getLab(c);
		for (short i=0; i<palette.length; ++i) {
			Color c2 = palette[i];			

			double curdist = sqr(c2.getAlpha() - c.getAlpha());
			if (curdist > mindist)
				continue;

			if (palette.length > 32) {
				curdist += PR * sqr(c2.getRed() - c.getRed());
				if (curdist > mindist)
					continue;

				curdist += PG * sqr(c2.getGreen() - c.getGreen());
				if (curdist > mindist)
					continue;

				curdist += PB * sqr(c2.getBlue() - c.getBlue());
			}
			else {
				Lab lab2 = getLab(c2);
				
				double deltaL_prime_div_k_L_S_L = CIELABConvertor.L_prime_div_k_L_S_L(lab1, lab2);
				curdist += sqr(deltaL_prime_div_k_L_S_L);
				if (curdist > mindist)
					continue;
	
				MutableDouble a1Prime = new MutableDouble(), a2Prime = new MutableDouble(), CPrime1 = new MutableDouble(), CPrime2 = new MutableDouble();
				double deltaC_prime_div_k_L_S_L = CIELABConvertor.C_prime_div_k_L_S_L(lab1, lab2, a1Prime, a2Prime, CPrime1, CPrime2);
				curdist += sqr(deltaC_prime_div_k_L_S_L);
				if (curdist > mindist)
					continue;
	
				MutableDouble barCPrime = new MutableDouble(), barhPrime = new MutableDouble();
				double deltaH_prime_div_k_L_S_L = CIELABConvertor.H_prime_div_k_L_S_L(lab1, lab2, a1Prime, a2Prime, CPrime1, CPrime2, barCPrime, barhPrime);
				curdist +=  sqr(deltaH_prime_div_k_L_S_L);
				if (curdist > mindist)
					continue;
	
				curdist += CIELABConvertor.R_T(barCPrime, barhPrime, deltaC_prime_div_k_L_S_L, deltaH_prime_div_k_L_S_L);
			}

			if (curdist > mindist)
				continue;
			mindist = curdist;
			k = i;
		}
		return k;
	}

	private short closestColorIndex(final Color[] palette, final Color c)
	{
		short k = 0;
		short[] closest = new short[5];
		short[] got = closestMap.get(c);
		if (got == null) {
			closest[2] = closest[3] = SHORT_MAX;
			Lab lab1 = getLab(c);

			for (; k < palette.length; k++) {
				Color c2 = palette[k];
				Lab lab2 = getLab(c2);

				closest[4] = (short) (sqr(lab2.alpha - lab1.alpha) + CIELABConvertor.CIEDE2000(lab2, lab1));
				//closest[4] = Math.abs(lab2.alpha - lab1.alpha) + Math.abs(lab2.L - lab1.L) + Math.abs(lab2.A - lab1.A) + Math.abs(lab2.B - lab1.B);
				if (closest[4] < closest[2]) {
					closest[1] = closest[0];
					closest[3] = closest[2];
					closest[0] = k;
					closest[2] = closest[4];
				}
				else if (closest[4] < closest[3]) {
					closest[1] = k;
					closest[3] = closest[4];
				}
			}

			if (closest[3] == SHORT_MAX)
				closest[2] = 0;
		}
		else
			closest = got;

		Random rand = new Random();
		if (closest[2] == 0 || (rand.nextInt(SHORT_MAX) % (closest[3] + closest[2])) <= closest[3])
			k = closest[0];
		else
			k = closest[1];

		closestMap.put(c, closest);
		return k;
	}

	@Override
	protected boolean quantize_image(final Color[] pixels, final Color[] palette, short[] qPixels, final boolean dither)
	{
		int nMaxColors = palette.length;

		int pixelIndex = 0;
		if (dither) {
			boolean odd_scanline = false;
			short[] row0, row1;
			int dir, k;
			final int DJ = 4;
			final int DITHER_MAX = 20;
			final int err_len = (width + 2) * DJ;
			int[] clamp = new int[DJ * 256];
			int[] limtb = new int[512];
			short[] erowerr = new short[err_len];
			short[] orowerr = new short[err_len];

			for (int i = 0; i < 256; i++) {
				clamp[i] = 0;
				clamp[i + 256] = (short) i;
				clamp[i + 512] = BYTE_MAX;
				clamp[i + 768] = BYTE_MAX;

				limtb[i] = -DITHER_MAX;
				limtb[i + 256] = DITHER_MAX;
			}
			for (int i = -DITHER_MAX; i <= DITHER_MAX; i++)
				limtb[i + 256] = i;

			for (short i = 0; i < height; i++) {
				if (odd_scanline) {
					dir = -1;
					pixelIndex += (width - 1);
					row0 = orowerr;
					row1 = erowerr;
				}
				else {
					dir = 1;
					row0 = erowerr;
					row1 = orowerr;
				}

				int cursor0 = DJ, cursor1 = width * DJ;
				row1[cursor1] = row1[cursor1 + 1] = row1[cursor1 + 2] = row1[cursor1 + 3] = 0;
				for (short j = 0; j < width; j++) {
					Color c = pixels[pixelIndex];
					int[] ditherPixel = calcDitherPixel(c, clamp, row0, cursor0, hasSemiTransparency);
					int r_pix = ditherPixel[0];
                    int g_pix = ditherPixel[1];
                    int b_pix = ditherPixel[2];
                    int a_pix = ditherPixel[3];

					Color c1 = new Color(r_pix, g_pix, b_pix, a_pix);
					qPixels[pixelIndex] = nearestColorIndex(palette, c1);

					Color c2 = palette[qPixels[pixelIndex]];
					if(nMaxColors > 256)
						qPixels[pixelIndex] = (short) getARGB1555(c2);

					r_pix = limtb[r_pix - c2.getRed() + 256];
					g_pix = limtb[g_pix - c2.getGreen() + 256];
					b_pix = limtb[b_pix - c2.getBlue() + 256];
					a_pix = limtb[a_pix - c2.getAlpha() + 256];

					k = r_pix * 2;
					row1[cursor1 - DJ] = (short) r_pix;
					row1[cursor1 + DJ] += (r_pix += k);
					row1[cursor1] += (r_pix += k);
					row0[cursor0 + DJ] += (r_pix += k);

					k = g_pix * 2;
					row1[cursor1 + 1 - DJ] = (short) g_pix;
					row1[cursor1 + 1 + DJ] += (g_pix += k);
					row1[cursor1 + 1] += (g_pix += k);
					row0[cursor0 + 1 + DJ] += (g_pix += k);

					k = b_pix * 2;
					row1[cursor1 + 2 - DJ] = (short) b_pix;
					row1[cursor1 + 2 + DJ] += (b_pix += k);
					row1[cursor1 + 2] += (b_pix += k);
					row0[cursor0 + 2 + DJ] += (b_pix += k);

					k = a_pix * 2;
					row1[cursor1 + 3 - DJ] = (short) a_pix;
					row1[cursor1 + 3 + DJ] += (a_pix += k);
					row1[cursor1 + 3] += (a_pix += k);
					row0[cursor0 + 3 + DJ] += (a_pix += k);

					cursor0 += DJ;
					cursor1 -= DJ;
					pixelIndex += dir;
				}
				if ((i % 2) == 1)
					pixelIndex += (width + 1);

				odd_scanline = !odd_scanline;
			}
			return true;
		}

		if(hasSemiTransparency || nMaxColors < 64) {
			for (int i = 0; i < qPixels.length; i++)
				qPixels[i] = nearestColorIndex(palette, pixels[i]);
		}
		else {
			for (int i = 0; i < qPixels.length; i++)
				qPixels[i] = closestColorIndex(palette, pixels[i]);
		}

		return true;
	}

	@Override
	public short[] convert(int nMaxColors, boolean dither) {
		final Color[] cPixels = new Color[pixels.length];		
		for (int i = 0; i<pixels.length; ++i) {
			int pixel = pixels[i];
			int alfa = (pixel >> 24) & 0xff;
			int r   = (pixel >> 16) & 0xff;
			int g = (pixel >>  8) & 0xff;
			int b  = (pixel      ) & 0xff;
			cPixels[i] = new Color(r, g, b, alfa);
			if (alfa < BYTE_MAX) {
				hasSemiTransparency = true;
				if (alfa == 0) {
					m_transparentPixelIndex = i;
					m_transparentColor = cPixels[i];
				}
			}			
		}

		if (hasSemiTransparency)
			PR = PG = PB = 1.0;
		Color[] palette = new Color[nMaxColors];
		boolean quan_sqrt = Math.random() < nMaxColors / 64.0;
		if (nMaxColors > 2)
			palette = pnnquan(cPixels, nMaxColors, quan_sqrt);
		else {
			if (hasSemiTransparency) {
				palette[0] = new Color(0, 0, 0, 0);
				palette[1] = Color.BLACK;
			}
			else {
				palette[0] = Color.BLACK;
				palette[1] = Color.WHITE;
			}
		}

		short[] qPixels = new short[cPixels.length];
		if (nMaxColors > 256)
			dither = true;
		quantize_image(cPixels, palette, qPixels, dither);
		if (m_transparentPixelIndex >= 0) {
			short k = qPixels[m_transparentPixelIndex];
			if (nMaxColors > 2)
				palette[k] = m_transparentColor;
			else if (!palette[k].equals(m_transparentColor)) {
				Color c1 = palette[0]; palette[0] = palette[1]; palette[1] = c1;
			}
		}
		pixelMap.clear();
		closestMap.clear();

		return qPixels;
	}

}
