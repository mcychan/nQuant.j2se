package nQuant.j2se;
/* Fast pairwise nearest neighbor based algorithm for multilevel thresholding
Copyright (C) 2004-2016 Mark Tyler and Dmitry Groshev
Copyright (c) 2018 Miller Cy Chan
 * error measure; time used is proportional to number of bins squared - WJ */

import java.awt.Color;
import java.awt.Image;
import java.awt.image.ColorModel;
import java.awt.image.DataBuffer;
import java.awt.image.DirectColorModel;
import java.awt.image.ImageObserver;
import java.awt.image.IndexColorModel;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

public class PnnQuantizer {
	protected final short SHORT_MAX = Short.MAX_VALUE;
	protected final char BYTE_MAX = -Byte.MIN_VALUE + Byte.MAX_VALUE;
	protected boolean hasSemiTransparency = false;
	protected int m_transparentPixelIndex = -1;
	protected final int width, height;	
	protected int pixels[] = null;
	protected Color m_transparentColor;
	protected ColorModel m_colorModel;
	protected Map<Color, short[]> closestMap = new HashMap<Color, short[]>();	

	public PnnQuantizer(Image im, int w, int h) throws IOException {
		width = w;
		height = h;
		setPixels(im);
	}

	public PnnQuantizer(Image im, ImageObserver obs) throws IOException {
		width = im.getWidth(obs);
		height = im.getHeight(obs);
		setPixels(im, obs);
	}

	private void setPixels(Image im, ImageObserver obs) throws IOException {
		if (im == null)
			throw new IOException ("Image is null");		
		setPixels(im);
	}

	private void setPixels(Image im) throws IOException {
		pixels = new int [width * height];
		java.awt.image.PixelGrabber pg
		= new java.awt.image.PixelGrabber(im, 0, 0, width, height, pixels, 0, width);
		try {
			pg.grabPixels();
		} catch (InterruptedException e) { }
		if ((pg.getStatus() & java.awt.image.ImageObserver.ABORT) != 0) {
			throw new IOException ("Image pixel grab aborted or errored");
		}
	}

	private static final class Pnnbin {
		double ac = 0, rc = 0, gc = 0, bc = 0, err = 0;
		int cnt = 0;
		int nn, fw, bk, tm, mtm;
	}

	protected int getColorIndex(final Color c, boolean hasSemiTransparency, int transparentPixelIndex )
	{
		if(hasSemiTransparency)
			return (c.getAlpha() & 0xF0) << 8 | (c.getRed() & 0xF0) << 4 | (c.getGreen() & 0xF0) | (c.getBlue() >> 4);
		if (transparentPixelIndex >= 0)
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
		int n1 = bin1.cnt;
		double wa = bin1.ac;
		double wr = bin1.rc;
		double wg = bin1.gc;
		double wb = bin1.bc;
		for (int i = bin1.fw; i != 0; i = bins[i].fw) {
			double nerr = sqr(bins[i].ac - wa) + sqr(bins[i].rc - wr) + sqr(bins[i].gc - wg) + sqr(bins[i].bc - wb);
			double n2 = bins[i].cnt;
			nerr *= (n1 * n2) / (n1 + n2);
			if (nerr >= err)
				continue;
			err = nerr;
			nn = i;
		}
		bin1.err = err;
		bin1.nn = nn;
	}

	protected void setColorModel(final List<Color> palette)
	{
		int nMaxColors = palette.size();

		int transparentARGB = -1;
		if(m_transparentColor != null)
			transparentARGB = (m_transparentColor.getAlpha() << 24) | (m_transparentColor.getRed() << 16) | (m_transparentColor.getGreen() << 8) | m_transparentColor.getBlue();
		
		if(nMaxColors <= 256) {
			int[] paletteARGB = new int[nMaxColors];
			for(int i=0; i<nMaxColors; ++i) {
				Color c1 = palette.get(i);
				paletteARGB[i] = (c1.getAlpha() << 24) | (c1.getRed() << 16) | (c1.getGreen() << 8) | c1.getBlue();
			}
			
			m_colorModel = new IndexColorModel(8,         // bits per pixel
				nMaxColors,         // size of color component array
				paletteARGB,   // color map
				0,         // offset in the map
				m_transparentPixelIndex >= 0,      // has alpha
				transparentARGB,         // the pixel value that should be transparent
				DataBuffer.TYPE_BYTE); // ARGB
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
	}

	private Color[] pnnquan(final Color[] pixels, int nMaxColors)
	{
		Pnnbin[] bins = new Pnnbin[65536];
		int[] heap = new int[65537];
		double err, n1, n2;

		/* Build histogram */
		for (final Color pixel : pixels) {
			// !!! Can throw gamma correction in here, but what to do about perceptual
			// !!! nonuniformity then?
			int index = getColorIndex(pixel, hasSemiTransparency, m_transparentPixelIndex);
			if(bins[index] == null)
				bins[index] = new Pnnbin();
			Pnnbin tb = bins[index];
			tb.ac += pixel.getAlpha();
			tb.rc += pixel.getRed();
			tb.gc += pixel.getGreen();
			tb.bc += pixel.getBlue();
			tb.cnt++;
		}

		/* Cluster nonempty bins at one end of array */
		int maxbins = 0;

		for (int i = 0; i < bins.length; ++i) {
			if (bins[i] == null)
				continue;

			double d = 1.0 / (double)bins[i].cnt;
			bins[i].ac *= d;
			bins[i].rc *= d;
			bins[i].gc *= d;
			bins[i].bc *= d;
			bins[i].cnt = (int) Math.sqrt(bins[i].cnt);
			bins[maxbins++] = bins[i];
		}

		for (int i = 0; i < maxbins - 1; i++) {
			bins[i].fw = (i + 1);
			bins[i + 1].bk = i;
		}
		// !!! Already zeroed out by calloc()
		//	bins[0].bk = bins[i].fw = 0;

		int h, l, l2 ;
		/* Initialize nearest neighbors and build heap of them */
		for (int i = 0; i < maxbins; i++) {
			find_nn(bins, i);
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
					find_nn(bins, b1);
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
			tb.rc = d * (n1 * tb.rc + n2 * nb.rc);
			tb.gc = d * (n1 * tb.gc + n2 * nb.gc);
			tb.bc = d * (n1 * tb.bc + n2 * nb.bc);
			tb.cnt += nb.cnt;
			tb.mtm = ++i;

			/* Unchain deleted bin */
			bins[nb.bk].fw = nb.fw;
			bins[nb.fw].bk = nb.bk;
			nb.mtm = 0xFFFF;
		}

		/* Fill palette */
		List<Color> palette = new ArrayList<Color>();
		short k = 0;
		for (int i = 0;; ++k) {
			int alpha = (int) Math.rint(bins[i].ac);
			palette.add(new Color((int) Math.rint(bins[i].rc), (int) Math.rint(bins[i].gc), (int) Math.rint(bins[i].bc), alpha));
			if (m_transparentPixelIndex >= 0 && palette.get(k).equals(m_transparentColor))
				Collections.swap(palette, 0, k);

			if ((i = bins[i].fw) == 0)
				break;
		}

		setColorModel(palette);
		return palette.toArray(new Color[0]);
	}

	private short nearestColorIndex(final Color[] palette, final int nMaxColors, final Color c)
	{
		short k = 0;
		double curdist, mindist = SHORT_MAX;
		for (int i=0; i<nMaxColors; ++i) {
			Color c2 = palette[i];

			double adist = Math.abs(c2.getAlpha() - c.getAlpha());
			curdist = adist;
			if (curdist > mindist)
				continue;

			double rdist = Math.abs(c2.getRed() - c.getRed());
			curdist += rdist;
			if (curdist > mindist)
				continue;

			double gdist = Math.abs(c2.getGreen() - c.getGreen());
			curdist += gdist;
			if (curdist > mindist)
				continue;

			double bdist = Math.abs(c2.getBlue() - c.getBlue());
			curdist += bdist;
			if (curdist > mindist)
				continue;

			mindist = curdist;
			k = (short) i;
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

			for (; k < palette.length; k++) {
				Color c2 = palette[k];

				closest[4] = (short) (Math.abs(c.getAlpha() - c2.getAlpha()) + Math.abs(c.getRed() - c2.getRed()) + Math.abs(c.getGreen() - c2.getGreen()) + Math.abs(c.getBlue() - c2.getBlue()));
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

	boolean quantize_image(final Color[] pixels, final Color[] palette, short[] qPixels, final boolean dither)
	{
		int nMaxColors = palette.length;

		int pixelIndex = 0;
		if (dither) {
			boolean odd_scanline = false;
			short[] row0, row1;
			int a_pix, r_pix, g_pix, b_pix, dir, k;
			final int DJ = 4;
			final int DITHER_MAX = 20;
			final int err_len = (width + 2) * DJ;
			int[] clamp = new int[DJ * 256];
			int[] limtb = new int[512];
			short[] erowerr = new short[err_len];
			short[] orowerr = new short[err_len];
			int[] lookup = new int[65536];

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
					r_pix = clamp[((row0[cursor0] + 0x1008) >> 4) + c.getRed()];
					g_pix = clamp[((row0[cursor0 + 1] + 0x1008) >> 4) + c.getGreen()];
					b_pix = clamp[((row0[cursor0 + 2] + 0x1008) >> 4) + c.getBlue()];
					a_pix = clamp[((row0[cursor0 + 3] + 0x1008) >> 4) + c.getAlpha()];

					Color c1 = new Color(r_pix, g_pix, b_pix, a_pix);
					int offset = getColorIndex(c1, hasSemiTransparency, m_transparentPixelIndex);
					if (lookup[offset] == 0)
						lookup[offset] = nearestColorIndex(palette, nMaxColors, c1) + 1;
					qPixels[pixelIndex] = (short) (lookup[offset] - 1);

					Color c2 = palette[qPixels[pixelIndex]];
					if(nMaxColors > 256)
						qPixels[pixelIndex] = (short) getColorIndex(c2, hasSemiTransparency, m_transparentPixelIndex);

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

		if(hasSemiTransparency || nMaxColors < 256) {
			for (int i = 0; i < qPixels.length; i++)
				qPixels[i] = nearestColorIndex(palette, nMaxColors, pixels[i]);
		}
		else {
			for (int i = 0; i < qPixels.length; i++)
				qPixels[i] = closestColorIndex(palette, pixels[i]);
		}

		return true;
	}

	public short[] convert (int nMaxColors, boolean dither) {
		final Color[] cPixels = new Color[pixels.length];		
		for (int i =0; i<pixels.length; ++i) {
			int pixel = pixels[i];
			int alfa = (pixel >> 24) & 0xff;
			int r   = (pixel >> 16) & 0xff;
			int g = (pixel >>  8) & 0xff;
			int b  = (pixel      ) & 0xff;
			cPixels[i] = new Color(pixels[i]);
			if (alfa < BYTE_MAX) {
				hasSemiTransparency = true;
				if (alfa == 0) {
					m_transparentPixelIndex = i;
					m_transparentColor = cPixels[i];
				}
			}			
		}

		Color[] palette = new Color[nMaxColors];
		if (nMaxColors > 2)
			palette = pnnquan(cPixels, nMaxColors);
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
		closestMap.clear();

		return qPixels;
	}

	public ColorModel getColorModel() {
		return m_colorModel;
	}

	public int getWidth() {
		return width;
	}

	public int getHeight() {
		return height;
	}		

}
