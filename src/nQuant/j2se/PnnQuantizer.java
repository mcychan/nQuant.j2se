package nQuant.j2se;
/* Fast pairwise nearest neighbor based algorithm for multilevel thresholding
Copyright (C) 2004-2016 Mark Tyler and Dmitry Groshev
Copyright (c) 2018-2020 Miller Cy Chan
 * error measure; time used is proportional to number of bins squared - WJ */

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.awt.image.ColorModel;
import java.awt.image.DataBuffer;
import java.awt.image.DirectColorModel;
import java.awt.image.ImageObserver;
import java.awt.image.IndexColorModel;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;

public class PnnQuantizer {
	protected final char BYTE_MAX = -Byte.MIN_VALUE + Byte.MAX_VALUE;
	protected short alphaThreshold = 0;
	protected boolean hasSemiTransparency = false;
	protected int m_transparentPixelIndex = -1;
	protected final int width, height;	
	protected int[] pixels;
	protected Color m_transparentColor;
	
	private double PR = .2126, PG = .7152, PB = .0722;
	private Color[] m_palette;
	protected ColorModel m_colorModel;
	protected Map<Integer, int[]> closestMap = new HashMap<>();
	protected Map<Integer, Short> nearestMap = new HashMap<>();

	private PnnQuantizer(BufferedImage im, int w, int h) throws IOException {
		width = w;
		height = h;
		setPixels(im);
	}

	public PnnQuantizer(BufferedImage im, ImageObserver obs) throws IOException {
		this(im, im.getWidth(obs), im.getHeight(obs));
	}

	private void setPixels(BufferedImage im) throws IOException {
	    pixels = im.getRGB(0, 0, width, height, null, 0, width);
	}

	private static final class Pnnbin {
		float ac = 0, rc = 0, gc = 0, bc = 0, err = 0;
		float cnt = 0;
		int nn, fw, bk, tm, mtm;
	}

	protected int getColorIndex(final Color c, boolean hasSemiTransparency, boolean hasTransparency)
	{
		if(hasSemiTransparency)
			return (c.getAlpha() & 0xF0) << 8 | (c.getRed() & 0xF0) << 4 | (c.getGreen() & 0xF0) | (c.getBlue() >> 4);
		if (hasTransparency)
			return (c.getAlpha() & 0x80) << 8 | (c.getRed() & 0xF8) << 7 | (c.getGreen() & 0xF8) << 2 | (c.getBlue() >> 3);
		return (c.getRed() & 0xF8) << 8 | (c.getGreen() & 0xFC) << 3 | (c.getBlue() >> 3);
	}
	
	protected int getColorIndex(final Color c)
	{
		if(hasSemiTransparency)
			return (c.getAlpha() & 0xF0) << 8 | (c.getRed() & 0xF0) << 4 | (c.getGreen() & 0xF0) | (c.getBlue() >> 4);
		if(m_transparentPixelIndex > -1)
			return (c.getAlpha() & 0x80) << 8 | (c.getRed() & 0xF8) << 7 | (c.getGreen() & 0xF8) << 2 | (c.getBlue() >> 3);
		return (c.getRed() & 0xF8) << 8 | (c.getGreen() & 0xFC) << 3 | (c.getBlue() >> 3);
	}

	protected double sqr(double value)
	{
		return value * value;
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
		for (int i = bin1.fw; i != 0; i = bins[i].fw) {
			double nerr = PR * sqr(bins[i].rc - wr) + PG * sqr(bins[i].gc - wg) + PB * sqr(bins[i].bc - wb);
			if(hasSemiTransparency)
				nerr += sqr(bins[i].ac - wa);
			
			float n2 = bins[i].cnt;
			nerr *= (n1 * n2) / (n1 + n2);
			if (nerr >= err)
				continue;
			err = nerr;
			nn = i;
		}
		bin1.err = (float) err;
		bin1.nn = nn;
	}
	
	private static int getBitsPerPixel(int nMaxColors)
	{
		return (int) Math.min(8, Math.log(nMaxColors) / Math.log(2));
	}
	
	protected final void setColorModel(final Color[] palette)
	{
		m_palette = palette;
		int nMaxColors = palette.length;
		
		if(nMaxColors <= 256) {
			int[] palettes = new int[nMaxColors];
			for(int i=0; i<nMaxColors; ++i) {
				Color c1 = palette[i];
				if(c1 == null)
					continue;

				palettes[i] = c1.getRGB();
			}
			
			m_colorModel = new IndexColorModel(getBitsPerPixel(nMaxColors),         // bits per pixel
				nMaxColors,         // size of color component array
				palettes,   // color map
                0,         // offset in the map
                m_transparentPixelIndex > -1,      // has alpha
                m_transparentPixelIndex,         // the pixel value that should be transparent
                DataBuffer.TYPE_BYTE);			
		}
		else if (hasSemiTransparency) {
			final int DCM_4444_RED_MASK = 0x0f00;
			final int DCM_4444_GRN_MASK = 0x00f0;
			final int DCM_4444_BLU_MASK = 0x000f;
			final int DCM_4444_ALP_MASK = 0xf000;
			m_colorModel = new DirectColorModel(16,
					DCM_4444_RED_MASK,
					DCM_4444_GRN_MASK,
					DCM_4444_BLU_MASK,
					DCM_4444_ALP_MASK);
		}
		else if (m_transparentPixelIndex >= 0) {
			final int DCM_1555_RED_MASK = 0x7c00;
			final int DCM_1555_GRN_MASK = 0x03e0;
			final int DCM_1555_BLU_MASK = 0x001f;
			final int DCM_1555_ALP_MASK = 0x8000;
			m_colorModel = new DirectColorModel(16,
					DCM_1555_RED_MASK,
					DCM_1555_GRN_MASK,
					DCM_1555_BLU_MASK,
					DCM_1555_ALP_MASK);
		}
		else {
			final int DCM_565_RED_MASK = 63488;
			final int DCM_565_GRN_MASK = 2016;
			final int DCM_565_BLU_MASK = 31;
			m_colorModel = new DirectColorModel(16,
					DCM_565_RED_MASK,
					DCM_565_GRN_MASK,
					DCM_565_BLU_MASK);
		}
	}

	protected Color[] pnnquan(final Color[] pixels, int nMaxColors, short quan_rt)
	{
		closestMap.clear();
		nearestMap.clear();
		Pnnbin[] bins = new Pnnbin[65536];		

		/* Build histogram */
		for (Color c : pixels) {
			// !!! Can throw gamma correction in here, but what to do about perceptual
			// !!! nonuniformity then?
			
			int index = getColorIndex(c, hasSemiTransparency, nMaxColors < 64 || m_transparentPixelIndex >= 0);
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
		if (sqr(nMaxColors) / maxbins < .03)
			quan_rt = 0;		

		int j = 0;
		for (; j < maxbins - 1; ++j) {
			bins[j].fw = j + 1;
			bins[j + 1].bk = j;
			
			if (quan_rt > 0)
				bins[j].cnt = (int) Math.sqrt(bins[j].cnt);				
			else if (quan_rt < 0)
				bins[j].cnt = (int) Math.cbrt(bins[j].cnt);
		}
		if (quan_rt > 0)
			bins[j].cnt = (int) Math.sqrt(bins[j].cnt);			
		else if (quan_rt < 0)
			bins[j].cnt = (int) Math.cbrt(bins[j].cnt);

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
			tb.ac = d * (float) Math.round(n1 * tb.ac + n2 * nb.ac);
			tb.rc = d * (float) Math.round(n1 * tb.rc + n2 * nb.rc);
			tb.gc = d * (float) Math.round(n1 * tb.gc + n2 * nb.gc);
			tb.bc = d * (float) Math.round(n1 * tb.bc + n2 * nb.bc);
			tb.cnt += nb.cnt;
			tb.mtm = ++i;

			/* Unchain deleted bin */
			bins[nb.bk].fw = nb.fw;
			bins[nb.fw].bk = nb.bk;
			nb.mtm = 0xFFFF;
		}

		/* Fill palette */
		Color[] palette = new Color[nMaxColors];
		short k = 0;
		for (int i = 0;; ++k) {
			int alpha = (int) bins[i].ac;
			palette[k] = new Color((int) bins[i].rc, (int) bins[i].gc, (int) bins[i].bc, alpha);
			if (m_transparentPixelIndex >= 0 && m_transparentColor.equals(palette[k])) {
				Color temp = palette[0]; palette[0] = palette[k]; palette[k] = temp;
			}

			if ((i = bins[i].fw) == 0)
				break;
		}

		setColorModel(palette);
		return palette;
	}

	public short nearestColorIndex(final Color[] palette, final Color c)
	{
		Short got = nearestMap.get(c.getRGB());
		if (got != null)
			return got;
		
		short k = 0;
		if (c.getAlpha() <= alphaThreshold)
            return k;
		
		double mindist = 1e100;
		for (int i = 0; i<palette.length; ++i) {
			Color c2 = palette[i];
			if(c2 == null)
				break;

			double curdist = Math.abs(c2.getAlpha() - c.getAlpha());
			if (curdist > mindist)
				continue;

			curdist += PR * Math.abs(c2.getRed() - c.getRed());
			if (curdist > mindist)
				continue;

			curdist += PG * Math.abs(c2.getGreen() - c.getGreen());
			if (curdist > mindist)
				continue;

			curdist += PB * Math.abs(c2.getBlue() - c.getBlue());
			if (curdist > mindist)
				continue;

			mindist = curdist;
			k = (short) i;
		}
		nearestMap.put(c.getRGB(), k);
		return k;
	}

	protected short closestColorIndex(final Color[] palette, final Color c)
	{
		short k = 0;
		int[] closest = new int[5];
		int[] got = closestMap.get(c.getRGB());
		if (got == null) {
			closest[2] = closest[3] = Integer.MAX_VALUE;

			for (; k < palette.length; ++k) {
				Color c2 = palette[k];
				if(c2 == null)
					break;

				closest[4] = (int) (Math.abs(c.getAlpha() - c2.getAlpha()) + Math.abs(c.getRed() - c2.getRed()) + Math.abs(c.getGreen() - c2.getGreen()) + Math.abs(c.getBlue() - c2.getBlue()));
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

			if (closest[3] == Integer.MAX_VALUE)
				closest[2] = 0;
		}
		else
			closest = got;

		Random rand = new Random();
		if (closest[2] == 0 || (rand.nextInt(32767) % (closest[3] + closest[2])) <= closest[3])
			k = (short) closest[0];
		else
			k = (short) closest[1];

		closestMap.put(c.getRGB(), closest);
		return k;
	}
	
	protected int[] calcDitherPixel(Color c, short[] clamp, int[] rowerr, int cursor, boolean noBias)
    {
        int[] ditherPixel = new int[4];
        if (noBias) {
            ditherPixel[0] = clamp[((rowerr[cursor] + 0x1008) >> 4) + c.getRed()];
            ditherPixel[1] = clamp[((rowerr[cursor + 1] + 0x1008) >> 4) + c.getGreen()];
            ditherPixel[2] = clamp[((rowerr[cursor + 2] + 0x1008) >> 4) + c.getBlue()];
            ditherPixel[3] = clamp[((rowerr[cursor + 3] + 0x1008) >> 4) + c.getAlpha()];
            return ditherPixel;
        }

        ditherPixel[0] = clamp[((rowerr[cursor] + 0x2010) >> 5) + c.getRed()];
        ditherPixel[1] = clamp[((rowerr[cursor + 1] + 0x1008) >> 4) + c.getGreen()];
        ditherPixel[2] = clamp[((rowerr[cursor + 2] + 0x2010) >> 5) + c.getBlue()];
        ditherPixel[3] = c.getAlpha();
        return ditherPixel;
    }
	
	protected short[] quantize_image(final Color[] pixels, final Color[] palette, final boolean dither)
	{
		short[] qPixels = new short[pixels.length];
		int nMaxColors = palette.length;

		int pixelIndex = 0;
		if (dither) {			
			final int DJ = 4;
			final int BLOCK_SIZE = 256;
			final short DITHER_MAX = 20;
			final int err_len = (width + 2) * DJ;
			short[] clamp = new short[DJ * BLOCK_SIZE];
			short[] limtb = new short[2 * BLOCK_SIZE];			

			for (short i = 0; i < BLOCK_SIZE; ++i) {
				clamp[i] = 0;
				clamp[i + BLOCK_SIZE] = i;
				clamp[i + BLOCK_SIZE * 2] = BYTE_MAX;
				clamp[i + BLOCK_SIZE * 3] = BYTE_MAX;

				limtb[i] = -DITHER_MAX;
				limtb[i + BLOCK_SIZE] = DITHER_MAX;
			}
			for (short i = -DITHER_MAX; i <= DITHER_MAX; ++i)
				limtb[i + BLOCK_SIZE] = i % 4 == 3 ? 0 : i;

			boolean noBias = hasSemiTransparency || nMaxColors < 64;
			int dir = 1;
			int[] row0 = new int[err_len];
			int[] row1 = new int[err_len];
			int[] lookup = new int[65536];
			for (int i = 0; i < height; ++i) {
				if (dir < 0)
					pixelIndex += width - 1;					

				int cursor0 = DJ, cursor1 = width * DJ;
				row1[cursor1] = row1[cursor1 + 1] = row1[cursor1 + 2] = row1[cursor1 + 3] = 0;
				for (int j = 0; j < width; ++j) {
					Color c = pixels[pixelIndex];
					int[] ditherPixel = calcDitherPixel(c, clamp, row0, cursor0, noBias);
					int r_pix = ditherPixel[0];
                    int g_pix = ditherPixel[1];
                    int b_pix = ditherPixel[2];
                    int a_pix = ditherPixel[3];

					Color c1 = new Color(r_pix, g_pix, b_pix, a_pix);
					if(noBias) {
						int offset = getColorIndex(c1, hasSemiTransparency, m_transparentPixelIndex >= 0);
						if (lookup[offset] == 0)
							lookup[offset] = (c.getAlpha() == 0) ? 1 : nearestColorIndex(palette, c1) + 1;
						qPixels[pixelIndex] = (short) (lookup[offset] - 1);
					}
					else
						qPixels[pixelIndex] = (c.getAlpha() == 0) ? 0 : nearestColorIndex(palette, c1);

					Color c2 = palette[qPixels[pixelIndex]];
					if(nMaxColors > 256)
						qPixels[pixelIndex] = (short) getColorIndex(c2);

					r_pix = limtb[r_pix - c2.getRed() + BLOCK_SIZE];
					g_pix = limtb[g_pix - c2.getGreen() + BLOCK_SIZE];
					b_pix = limtb[b_pix - c2.getBlue() + BLOCK_SIZE];
					a_pix = limtb[a_pix - c2.getAlpha() + BLOCK_SIZE];

					int k = r_pix * 2;
					row1[cursor1 - DJ] = r_pix;
					row1[cursor1 + DJ] += (r_pix += k);
					row1[cursor1] += (r_pix += k);
					row0[cursor0 + DJ] += (r_pix + k);

					k = g_pix * 2;
					row1[cursor1 + 1 - DJ] = g_pix;
					row1[cursor1 + 1 + DJ] += (g_pix += k);
					row1[cursor1 + 1] += (g_pix += k);
					row0[cursor0 + 1 + DJ] += (g_pix + k);

					k = b_pix * 2;
					row1[cursor1 + 2 - DJ] = b_pix;
					row1[cursor1 + 2 + DJ] += (b_pix += k);
					row1[cursor1 + 2] += (b_pix += k);
					row0[cursor0 + 2 + DJ] += (b_pix + k);

					k = a_pix * 2;
					row1[cursor1 + 3 - DJ] = a_pix;
					row1[cursor1 + 3 + DJ] += (a_pix += k);
					row1[cursor1 + 3] += (a_pix += k);
					row0[cursor0 + 3 + DJ] += (a_pix + k);

					cursor0 += DJ;
					cursor1 -= DJ;
					pixelIndex += dir;
				}
				if ((i % 2) == 1)
					pixelIndex += width + 1;

				dir *= -1;
				int[] temp = row0; row0 = row1; row1 = temp;
			}
			return qPixels;
		}

		if(hasSemiTransparency || nMaxColors < 64) {
			for (int i = 0; i < qPixels.length; ++i)
				qPixels[i] = nearestColorIndex(palette, pixels[i]);
		}
		else {
			for (int i = 0; i < qPixels.length; ++i)
				qPixels[i] = closestColorIndex(palette, pixels[i]);
		}

		return qPixels;
	}
	
	protected Ditherable getDitherFn(final boolean dither) {
		return new Ditherable() {
			@Override
			public int getColorIndex(Color c) {
				return PnnQuantizer.this.getColorIndex(c, hasSemiTransparency, m_transparentPixelIndex >= 0);
			}
			
			@Override
			public short nearestColorIndex(Color[] palette, Color c) {
				if(dither)
					return PnnQuantizer.this.nearestColorIndex(palette, c);
				return PnnQuantizer.this.closestColorIndex(palette, c);
			}
			
		};
	}
	
	protected short[] dither(final Color[] cPixels, Color[] palette, int nMaxColors, int width, int height, boolean dither)
    {
		short[] qPixels;
		if ((nMaxColors < 64 && nMaxColors > 32) || hasSemiTransparency)
			qPixels = quantize_image(cPixels, palette, dither);
		else if(nMaxColors <= 32)
			qPixels = GilbertCurve.dither(width, height, cPixels, palette, getDitherFn(dither), nMaxColors > 2 ? 1.8f : 1.5f);
		else			
			qPixels = GilbertCurve.dither(width, height, cPixels, palette, getDitherFn(dither));	
		
		if (!dither)
			BlueNoise.dither(width, height, cPixels, palette, getDitherFn(dither), qPixels, 1.0f);
			
		closestMap.clear();
		nearestMap.clear();
		return qPixels;
    }

	public short[] convert(int nMaxColors, boolean dither) {
		final Color[] cPixels = new Color[pixels.length];		
		for (int i = 0; i<pixels.length; ++i) {
			int pixel = pixels[i];
			int alfa = (pixel >> 24) & 0xff;
			cPixels[i] = new Color(pixel, true);
			if (alfa < BYTE_MAX) {				
				if (alfa == 0) {
					m_transparentPixelIndex = i;
					m_transparentColor = cPixels[i];
				}
				else
					hasSemiTransparency = true;
			}			
		}

		if (hasSemiTransparency || nMaxColors <= 32)
            PR = PG = PB = 1;
		else if(width < 512 || height < 512) {
			PR = 0.299; PG = 0.587; PB = 0.114;
		}
		
		Color[] palette;
		if (nMaxColors > 2)
			palette = pnnquan(cPixels, nMaxColors, (short) 1);
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
			setColorModel(palette);	
		}

		if (nMaxColors > 256)
			dither = true;

		short[] qPixels = dither(cPixels, palette, nMaxColors, width, height, dither);

		if (m_transparentPixelIndex >= 0) {
			short k = qPixels[m_transparentPixelIndex];
			if (nMaxColors > 2)
				palette[k] = m_transparentColor;
			else if (!palette[k].equals(m_transparentColor)) {
				Color c1 = palette[0]; palette[0] = palette[1]; palette[1] = c1;
			}
		}		

		return qPixels;
	}

	public ColorModel getColorModel() {
		return m_colorModel;
	}	

	public Color[] getPalette() {
		return m_palette;
	}

	public int getWidth() {
		return width;
	}

	public int getHeight() {
		return height;
	}
	
	public boolean hasAlpha() {
		return m_transparentPixelIndex > -1;
	}

}
