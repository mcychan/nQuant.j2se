package nQuant.j2se;
/* Fast pairwise nearest neighbor based algorithm for multilevel thresholding
Copyright (C) 2004-2016 Mark Tyler and Dmitry Groshev
Copyright (c) 2018-2025 Miller Cy Chan
 * error measure; time used is proportional to number of bins squared - WJ */

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.awt.image.ColorModel;
import java.awt.image.ImageObserver;
import java.util.HashMap;
import java.util.Map;

import static nQuant.j2se.BitmapUtilities.BYTE_MAX;

public class PnnQuantizer {
	protected short alphaThreshold = 0xF;
	protected boolean hasSemiTransparency = false;
	protected int m_transparentPixelIndex = -1;
	protected final int width, height;	
	protected int[] pixels;
	protected Color m_transparentColor = new Color(BYTE_MAX, BYTE_MAX, BYTE_MAX, 0);
	
	protected double PR = 0.299, PG = 0.587, PB = 0.114, PA = .3333;	
	protected double ratio = .5, weight = 1;
	protected static final float[][] coeffs = new float[][] {
		{0.299f, 0.587f, 0.114f},
		{-0.14713f, -0.28886f, 0.436f},
		{0.615f, -0.51499f, -0.10001f}
	};
	
	protected Color[] m_palette;
	private ColorModel m_colorModel;
	protected Map<Integer, int[]> closestMap = new HashMap<>();
	protected Map<Integer, Short> nearestMap = new HashMap<>();	

	private PnnQuantizer(BufferedImage im, int w, int h) {
		width = w;
		height = h;
		setPixels(im);
	}
	
	protected PnnQuantizer(PnnQuantizer quantizer) {
		alphaThreshold = quantizer.alphaThreshold;
		hasSemiTransparency = quantizer.hasSemiTransparency;
		m_transparentPixelIndex = quantizer.m_transparentPixelIndex;
		width = quantizer.width;
		height = quantizer.height;
		pixels = quantizer.pixels;
		m_transparentColor = quantizer.m_transparentColor;
	}

	public PnnQuantizer(BufferedImage im, ImageObserver obs) {
		this(im, im.getWidth(obs), im.getHeight(obs));
	}

	private void setPixels(BufferedImage im) {
		pixels = im.getRGB(0, 0, width, height, null, 0, width);
	}

	private static final class Pnnbin {
		double ac = 0, rc = 0, gc = 0, bc = 0;
		float cnt = 0, err = 0;
		int nn, fw, bk, tm, mtm;
	}	
	
	protected int getColorIndex(final Color c)
	{
		if(hasSemiTransparency)
			return (c.getAlpha() & 0xF0) << 8 | (c.getRed() & 0xF0) << 4 | (c.getGreen() & 0xF0) | (c.getBlue() >> 4);
		if(m_transparentPixelIndex > -1)
			return (c.getAlpha() & 0x80) << 8 | (c.getRed() & 0xF8) << 7 | (c.getGreen() & 0xF8) << 2 | (c.getBlue() >> 3);
		return (c.getRed() & 0xF8) << 8 | (c.getGreen() & 0xFC) << 3 | (c.getBlue() >> 3);
	}
	

	private void find_nn(Pnnbin[] bins, int idx)
	{
		int nn = 0;
		double err = 1e100;

		Pnnbin bin1 = bins[idx];
		float n1 = bin1.cnt;
		double wa = bin1.ac;
		double wr = bin1.rc;
		double wg = bin1.gc;
		double wb = bin1.bc;
		
		int start = 0;
		if(BlueNoise.TELL_BLUE_NOISE[idx & 4095] > 0)
			start = (PG < coeffs[0][1]) ? coeffs.length : 1;
		
		for (int i = bin1.fw; i != 0; i = bins[i].fw) {
			double n2 = bins[i].cnt, nerr2 = (n1 * n2) / (n1 + n2);
			if (nerr2 >= err)
				continue;
			
			double nerr = 0.0;
			if(hasSemiTransparency) {
				nerr += nerr2 * PA * BitmapUtilities.sqr(bins[i].ac - wa);
				if (nerr >= err)
					continue;
			}
			
			nerr += nerr2 * (1 - ratio) * PR * BitmapUtilities.sqr(bins[i].rc - wr);
			if (nerr >= err)
				continue;

			nerr += nerr2 * (1 - ratio) * PG * BitmapUtilities.sqr(bins[i].gc - wg);
			if (nerr >= err)
				continue;

			nerr += nerr2 * (1 - ratio) * PB * BitmapUtilities.sqr(bins[i].bc - wb);
			if (nerr >= err)
				continue;
			
			for (int j = start; j < coeffs.length; ++j) {
				nerr += nerr2 * ratio * BitmapUtilities.sqr(coeffs[j][0] * (bins[i].rc - wr));
				if (nerr >= err)
					break;
				
				nerr += nerr2 * ratio * BitmapUtilities.sqr(coeffs[j][1] * (bins[i].gc - wg));
				if (nerr >= err)
					break;
				
				nerr += nerr2 * ratio * BitmapUtilities.sqr(coeffs[j][2] * (bins[i].bc - wb));
				if (nerr >= err)
					break;
			}
			
			err = nerr;
			nn = i;
		}
		bin1.err = (float) err;
		bin1.nn = nn;
	}
	
	protected final void setColorModel(final Color[] palette)
	{
		m_palette = palette;
		m_colorModel = BitmapUtilities.createColorModel(palette, m_transparentPixelIndex, hasSemiTransparency, hasSemiTransparency);
	}

	@FunctionalInterface
	protected interface QuanFn {
		float get(float cnt);
	}

	protected QuanFn getQuanFn(int nMaxColors, short quan_rt) {
		if (quan_rt > 0) {
			if (nMaxColors < 64)
				return cnt -> (float) Math.sqrt(cnt);
			return cnt -> (int) Math.sqrt(cnt);
		}
		if (quan_rt < 0)
			return cnt -> (int) Math.cbrt(cnt);
		return cnt -> cnt;
	}

	protected Color[] pnnquan(final Color[] pixels, int nMaxColors)
	{
		closestMap.clear();
		nearestMap.clear();
		short quan_rt = (short) 1;
		Pnnbin[] bins = new Pnnbin[65536];

		/* Build histogram */
		for (int i = 0; i < pixels.length; ++i) {
			Color c = pixels[i];
			if (c.getAlpha() <= alphaThreshold)
				c = m_transparentColor;
			
			int index = BitmapUtilities.getColorIndex(c, hasSemiTransparency, nMaxColors < 64 || hasAlpha());
			if(bins[index] == null)
				bins[index] = new Pnnbin();
			Pnnbin tb = bins[index];			
			tb.ac += c.getAlpha();
			tb.rc += c.getRed();
			tb.gc += c.getGreen();
			tb.bc += c.getBlue();
			tb.cnt += 1.0f;
		}

		/* Cluster nonempty bins at one end of array */
		int maxbins = 0;
		for (int i = 0; i < bins.length; ++i) {
			if (bins[i] == null)
				continue;

			float d = 1f / bins[i].cnt;
			bins[i].ac *= d;
			bins[i].rc *= d;
			bins[i].gc *= d;
			bins[i].bc *= d;

			bins[maxbins++] = bins[i];
		}
		
		if (nMaxColors < 16)
			quan_rt = -1;
		
		weight = nMaxColors * 1.0 / maxbins;
		if (weight > .003 && weight < .005)
			quan_rt = 0;
		if (weight < .04 && PG >= coeffs[0][1]) {
			PR = PG = PB = PA = 1;
			if (nMaxColors >= 64)
				quan_rt = 0;
		}
		
		QuanFn quanFn = getQuanFn(nMaxColors, quan_rt);

		int j = 0;
		for (; j < maxbins - 1; ++j) {
			bins[j].fw = j + 1;
			bins[j + 1].bk = j;
			
			bins[j].cnt = quanFn.get(bins[j].cnt);
		}
		bins[j].cnt = quanFn.get(bins[j].cnt);

		int h, l, l2 ;
		/* Initialize nearest neighbors and build heap of them */
		int[] heap = new int[bins.length + 1];
		for (int i = 0; i < maxbins; ++i) {
			find_nn(bins, i);
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
					find_nn(bins, b1);
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
			float d = 1f / (n1 + n2);
			tb.ac = d * Math.round(n1 * tb.ac + n2 * nb.ac);
			tb.rc = d * Math.round(n1 * tb.rc + n2 * nb.rc);
			tb.gc = d * Math.round(n1 * tb.gc + n2 * nb.gc);
			tb.bc = d * Math.round(n1 * tb.bc + n2 * nb.bc);
			tb.cnt += n2;
			tb.mtm = ++i;

			/* Unchain deleted bin */
			bins[nb.bk].fw = nb.fw;
			bins[nb.fw].bk = nb.bk;
			nb.mtm = 0xFFFF;
		}

		/* Fill palette */
		Color[] palette = new Color[extbins > 0 ? nMaxColors : maxbins];
		short k = 0;
		for (int i = 0; k < palette.length; ++k) {
			palette[k] = new Color((int) bins[i].rc, (int) bins[i].gc, (int) bins[i].bc, (int) bins[i].ac);			

			i = bins[i].fw;
		}

		return palette;
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
		
		double pr = PR, pg = PG, pb = PB;
		
		double mindist = 1e100;
		for (int i = k; i < palette.length; ++i) {
			Color c2 = palette[i];

			double curdist = PA * BitmapUtilities.sqr(c2.getAlpha() - c.getAlpha());
			if (curdist > mindist)
				continue;

			curdist += pr * BitmapUtilities.sqr(c2.getRed() - c.getRed());
			if (curdist > mindist)
				continue;

			curdist += pg * BitmapUtilities.sqr(c2.getGreen() - c.getGreen());
			if (curdist > mindist)
				continue;

			curdist += pb * BitmapUtilities.sqr(c2.getBlue() - c.getBlue());
			if (curdist > mindist)
				continue;

			mindist = curdist;
			k = (short) i;
		}
		nearestMap.put(c.getRGB(), k);
		return k;
	}

	protected short closestColorIndex(final Color[] palette, Color c, final int pos)
	{
		if (c.getAlpha() <= alphaThreshold)
			return nearestColorIndex(palette, c, pos);
		
		int[] closest = new int[4];
		int[] got = closestMap.get(c.getRGB());
		if (got == null) {
			closest[2] = closest[3] = Integer.MAX_VALUE;
			
			double pr = PR, pg = PG, pb = PB;

			for (short k = 0; k < palette.length; ++k) {
				Color c2 = palette[k];

				double err = pr * BitmapUtilities.sqr(c2.getRed() - c.getRed());
				if (err >= closest[3])
					continue;
				
				err += pg * BitmapUtilities.sqr(c2.getGreen() - c.getGreen());
				if (err >= closest[3])
					continue;
				
				err += pb * BitmapUtilities.sqr(c2.getBlue() - c.getBlue());
				if (err >= closest[3])
					continue;
				
				if(hasSemiTransparency)
					err += PA * BitmapUtilities.sqr(c2.getAlpha() - c.getAlpha());
				
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
			
			closestMap.put(c.getRGB(), closest);
		}
		else
			closest = got;

		int MAX_ERR = palette.length << 2;
		int idx = (pos + 1) % 2;
		if (closest[3] * .67 < (closest[3] - closest[2]))
			idx = 0;
		else if (closest[0] > closest[1])
			idx = pos % 2;
			
		if(closest[idx + 2] >= MAX_ERR || (hasAlpha() && closest[idx + 2] == 0))
			return nearestColorIndex(palette, c, pos);
		return (short) closest[idx];
	}		
	
	protected Ditherable getDitherFn(final boolean dither) {
		return new Ditherable() {
			@Override
			public int getColorIndex(Color c) {
				return BitmapUtilities.getColorIndex(c, hasSemiTransparency, hasAlpha());
			}
			
			@Override
			public short nearestColorIndex(Color[] palette, Color c, final int pos) {
				if(dither)
					return PnnQuantizer.this.nearestColorIndex(palette, c, pos);
				return PnnQuantizer.this.closestColorIndex(palette, c, pos);
			}
			
		};
	}
	
	protected short[] dither(Color[] palette, int width, int height, boolean dither)
	{
		Ditherable ditherable = getDitherFn(dither);
		if(hasSemiTransparency)
			weight *= -1;
		short[] qPixels = GilbertCurve.dither(width, height, pixels, palette, ditherable, null, weight, dither);

		if (!dither && palette.length > 32)
			BlueNoise.dither(width, height, pixels, palette, ditherable, qPixels, 1.0f);
		
		closestMap.clear();
		nearestMap.clear();
		return qPixels;
	}

	protected Color[] grabPixels(int[] pixels, int nMaxColors, ThreadLocal<Boolean> isSemiTransparency) {
		final Color[] cPixels = new Color[pixels.length];
		int semiTransCount = 0;
		for (int i = 0; i < pixels.length; ++i) {
			int pixel = pixels[i];
			int alfa = (pixel >> 24) & 0xff;
			cPixels[i] = new Color(pixel, true);
			if(m_transparentPixelIndex > -1 && alfa == 0)
				cPixels[i] = m_transparentColor;
			
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

	public BufferedImage convert(int nMaxColors, boolean dither) {
		if (nMaxColors <= 32)
			PR = PG = PB = PA = 1;
		else {
			PR = coeffs[0][0]; PG = coeffs[0][1]; PB = coeffs[0][2];
		}
		
		Color[] palette;
		if(m_palette == null) {
			final Color[] cPixels = grabPixels(pixels, nMaxColors, null);
			if (nMaxColors > 2)
				palette = pnnquan(cPixels, nMaxColors);
			else {
				palette = new Color[nMaxColors];
				if (hasAlpha()) {
					palette[0] = m_transparentColor;
					palette[1] = Color.BLACK;
				}
				else {
					palette[0] = Color.BLACK;
					palette[1] = Color.WHITE;
				}
			}
		}
		else
			palette = m_palette;
			
		short[] qPixels = dither(palette, width, height, dither);
		if (hasAlpha() && nMaxColors > 2)
		{
			short k = qPixels[m_transparentPixelIndex];
			palette[k] = m_transparentColor;
		}
		setColorModel(palette);
		return BitmapUtilities.processImagePixels(qPixels, m_colorModel, width, height);
	}

	public Color[] getPalette() {
		return m_palette;
	}
	
	public boolean hasAlpha() {
		return m_transparentPixelIndex > -1;
	}

}
