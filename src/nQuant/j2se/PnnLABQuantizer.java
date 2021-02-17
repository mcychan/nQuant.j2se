package nQuant.j2se;

import java.awt.Color;
import java.awt.Image;
import java.awt.image.ImageObserver;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;

import nQuant.j2se.CIELABConvertor.Lab;
import nQuant.j2se.CIELABConvertor.MutableDouble;

public class PnnLABQuantizer extends PnnQuantizer {
	private double PR = .299, PG = .587, PB = .114;
	private double ratio = 1.0;
	private Map<Integer, Lab> pixelMap = new HashMap<Integer, Lab>();	

	public PnnLABQuantizer(Image im, int w, int h) throws IOException {
		super(im, w, h);
	}

	public PnnLABQuantizer(Image im, ImageObserver obs) throws IOException {
		super(im, obs);
	}

	private static final class Pnnbin {
		float ac = 0, Lc = 0, Ac = 0, Bc = 0, err = 0;
		int cnt = 0;
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
		int n1 = bin1.cnt;
		
		Lab lab1 = new Lab();
		lab1.alpha = bin1.ac; lab1.L = bin1.Lc; lab1.A = bin1.Ac; lab1.B = bin1.Bc;
		for (int i = bin1.fw; i != 0; i = bins[i].fw) {
			float n2 = bins[i].cnt;
			double nerr2 = (n1 * n2) / (n1 + n2);
			if (nerr2 >= err)
				continue;
			
			Lab lab2 = new Lab();
			lab2.alpha = bins[i].ac; lab2.L = bins[i].Lc; lab2.A = bins[i].Ac; lab2.B = bins[i].Bc;
			float alphaDiff = hasSemiTransparency ? Math.abs(lab2.alpha - lab1.alpha) : 0;
			double nerr = nerr2 * sqr(alphaDiff) * alphaDiff / 3.0;
			if (nerr >= err)
				continue;
			
			nerr += (1 - ratio) * nerr2 * sqr(lab2.L - lab1.L);
			if (nerr >= err)
				continue;

			nerr += (1 - ratio) * nerr2 * sqr(lab2.A - lab1.A);
			if (nerr >= err)
				continue;

			nerr += (1 - ratio) * nerr2 * sqr(lab2.B - lab1.B);

			if (nerr > err)
				continue;
			
			float deltaL_prime_div_k_L_S_L = CIELABConvertor.L_prime_div_k_L_S_L(lab1, lab2);
			nerr += ratio * nerr2 * sqr(deltaL_prime_div_k_L_S_L);
			if (nerr > err)
				continue;

			MutableDouble a1Prime = new MutableDouble(), a2Prime = new MutableDouble(), CPrime1 = new MutableDouble(), CPrime2 = new MutableDouble();
			float deltaC_prime_div_k_L_S_L = CIELABConvertor.C_prime_div_k_L_S_L(lab1, lab2, a1Prime, a2Prime, CPrime1, CPrime2);
			nerr += ratio * nerr2 * sqr(deltaC_prime_div_k_L_S_L);
			if (nerr > err)
				continue;

			MutableDouble barCPrime = new MutableDouble(), barhPrime = new MutableDouble();
			float deltaH_prime_div_k_L_S_L = CIELABConvertor.H_prime_div_k_L_S_L(lab1, lab2, a1Prime, a2Prime, CPrime1, CPrime2, barCPrime, barhPrime);
			nerr += ratio * nerr2 * sqr(deltaH_prime_div_k_L_S_L);
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

	private Color[] pnnquan(final Color[] pixels, int nMaxColors, boolean quan_sqrt)
	{
		Pnnbin[] bins = new Pnnbin[65536];
		int[] heap = new int[65537];

		/* Build histogram */
		for (final Color c : pixels) {
			// !!! Can throw gamma correction in here, but what to do about perceptual
			// !!! nonuniformity then?
			int index = getColorIndex(c, hasSemiTransparency);
			Lab lab1 = getLab(c.getRGB());
			if(bins[index] == null)
				bins[index] = new Pnnbin();
			Pnnbin tb = bins[index];
			tb.ac += c.getAlpha();
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

			float d = 1f / (float)bins[i].cnt;
			bins[i].ac *= d;
			bins[i].Lc *= d;
			bins[i].Ac *= d;
			bins[i].Bc *= d;
			
			if(quan_sqrt)
				bins[i].cnt = (int) Math.pow(bins[i].cnt, 0.6);
			bins[maxbins++] = bins[i];
		}

		for (int i = 0; i < maxbins - 1; i++) {
			bins[i].fw = (i + 1);
			bins[i + 1].bk = i;
		}
		// !!! Already zeroed out by calloc()
		//	bins[0].bk = bins[i].fw = 0;

		int h, l, l2;
		if(nMaxColors < 64)
			ratio = Math.min(1.0, Math.pow(nMaxColors, 2.1) / maxbins);
		else
			ratio = Math.min(1.0, Math.pow(nMaxColors, 1.05) / pixelMap.size());
		/* Initialize nearest neighbors and build heap of them */
		for (int i = 0; i < maxbins; i++) {
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
						l2++;
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
		for (int i = 0;; ++k) {
			Lab lab1 = new Lab();
			lab1.alpha = (int) bins[i].ac;
			lab1.L = bins[i].Lc; lab1.A = bins[i].Ac; lab1.B = bins[i].Bc;
			palette[k] = CIELABConvertor.LAB2RGB(lab1);
			if (m_transparentPixelIndex >= 0 && m_transparentColor.equals(palette[k])) {
				Color temp = palette[0]; palette[0] = palette[k]; palette[k] = temp;
			}

			if ((i = bins[i].fw) == 0)
				break;
		}

		setColorModel(palette);		
		return palette;
	}

	@Override
	protected short nearestColorIndex(final Color[] palette, final Color c)
	{
		Short got = nearestMap.get(c.getRGB());
		if (got != null)
			return got;
		
		short k = 0;
		double mindist = SHORT_MAX;
		Lab lab1 = getLab(c.getRGB());
		for (short i=0; i<palette.length; ++i) {
			Color c2 = palette[i];
			if(c2 == null)
				break;

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
				if (PB < 1) {
					if (curdist > mindist)
						continue;

					double luma1 = c.getRed() * PR + c.getGreen() * PG + c.getBlue() * PB;
					double luma2 = c2.getRed() * PR + c2.getGreen() * PG + c2.getBlue() * PB;
                    curdist += sqr(luma1 - luma2) / 3.0;
				}
			}
			else {		
				Lab lab2 = getLab(c2.getRGB());
				double deltaL_prime_div_k_L_S_L = CIELABConvertor.L_prime_div_k_L_S_L(lab1, lab2);
				curdist += sqr(deltaL_prime_div_k_L_S_L);
				if (curdist > mindist)
					continue;
	
				MutableDouble a1Prime = new MutableDouble(), a2Prime = new MutableDouble(), CPrime1 = new MutableDouble(), CPrime2 = new MutableDouble();
				float deltaC_prime_div_k_L_S_L = CIELABConvertor.C_prime_div_k_L_S_L(lab1, lab2, a1Prime, a2Prime, CPrime1, CPrime2);
				curdist += sqr(deltaC_prime_div_k_L_S_L);
				if (curdist > mindist)
					continue;
	
				MutableDouble barCPrime = new MutableDouble(), barhPrime = new MutableDouble();
				float deltaH_prime_div_k_L_S_L = CIELABConvertor.H_prime_div_k_L_S_L(lab1, lab2, a1Prime, a2Prime, CPrime1, CPrime2, barCPrime, barhPrime);
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
		nearestMap.put(c.getRGB(), k);
		return k;
	}

	@Override
	protected short closestColorIndex(final Color[] palette, final Color c)
	{
		short k = 0;
		short[] closest = new short[5];
		short[] got = closestMap.get(c.getRGB());
		if (got == null) {
			closest[2] = closest[3] = SHORT_MAX;
			Lab lab1 = getLab(c.getRGB());

			for (; k < palette.length; k++) {
				Color c2 = palette[k];
				Lab lab2 = getLab(c2.getRGB());

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

		closestMap.put(c.getRGB(), closest);
		return k;
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
		Color[] palette;
		boolean quan_sqrt = true;
		if (nMaxColors > 2)
			palette = pnnquan(cPixels, nMaxColors, quan_sqrt);
		else {
			palette = new Color[nMaxColors];
			if (hasSemiTransparency) {
				palette[0] = new Color(0, 0, 0, 0);
				palette[1] = Color.BLACK;
			}
			else {
				palette[0] = Color.BLACK;
				palette[1] = Color.WHITE;
			}
		}

		if (nMaxColors > 256)
			dither = true;
		short[] qPixels = quantize_image(cPixels, palette, dither);
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
		nearestMap.clear();

		return qPixels;
	}

}