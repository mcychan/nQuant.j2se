package nQuant.j2se;
/* Fast pairwise nearest neighbor based algorithm for multilevel thresholding
Copyright (C) 2004-2016 Mark Tyler and Dmitry Groshev
Copyright (c) 2018 Miller Cy Chan
* error measure; time used is proportional to number of bins squared - WJ */

import java.awt.Color;
import java.awt.Image;
import java.awt.image.ImageObserver;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class PnnQuant {
	private final short SHORT_MAX = Short.MAX_VALUE;
	private final char BYTE_MAX = -Byte.MIN_VALUE + Byte.MAX_VALUE;
	private boolean hasTransparency = false, hasSemiTransparency = false;
	private int pixels[] = null;
	private Color m_transparentColor;
	private Map<Color, short[]> closestMap = new HashMap<>();	
	
	public PnnQuant(Image im, int w, int h) throws IOException {
		setPixels(im, w, h);
	}

	public PnnQuant(Image im, ImageObserver obs) throws IOException {
		setPixels(im, obs);
	}
	
	private void setPixels(Image im, ImageObserver obs) throws IOException {
		if (im == null)
			throw new IOException ("Image is null");
		int w = im.getWidth(obs);
		int h = im.getHeight(obs);
		setPixels(im, w, h);
	}

	private void setPixels(Image im, int w, int h) throws IOException {
		pixels = new int [w * h];
		java.awt.image.PixelGrabber pg
		= new java.awt.image.PixelGrabber(im, 0, 0, w, h, pixels, 0, w);
		try {
			pg.grabPixels();
		} catch (InterruptedException e) { }
		if ((pg.getStatus() & java.awt.image.ImageObserver.ABORT) != 0) {
			throw new IOException ("Image pixel grab aborted or errored");
		}
	}

	private static class Pnnbin {
		double ac = 0, rc = 0, gc = 0, bc = 0, err = 0;
		int cnt = 0;
		int nn, fw, bk, tm, mtm;
	}

	private int getColorIndex(final Color c, boolean semiTransparency)
	{
		if(semiTransparency)
			return (c.getAlpha() & 0xF0) << 8 | (c.getRed() & 0xF0) << 4 | (c.getGreen() & 0xF0) | (c.getBlue() >> 4);
		if (hasTransparency)
			return (c.getAlpha() & 0x80) << 8 | (c.getRed() & 0xF8) << 7 | (c.getGreen() & 0xF8) << 2 | (c.getBlue() >> 3);
		return (c.getRed() & 0xF8) << 8 | (c.getGreen() & 0xFC) << 3 | (c.getBlue() >> 3);
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
			double nerr = Math.pow((bins[i].ac - wa), 2) + Math.pow((bins[i].rc - wr), 2) + Math.pow((bins[i].gc - wg), 2) + Math.pow((bins[i].bc - wb), 2);
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

	private Color[] pnnquan(final Color[] pixels, Pnnbin[] bins, int nMaxColors, boolean quan_sqrt)
	{
		int[] heap = new int[65537];
		double err, n1, n2;
		int l, l2, h, b1, maxbins, extbins;

		/* Build histogram */
		for (final Color pixel : pixels) {
			// !!! Can throw gamma correction in here, but what to do about perceptual
			// !!! nonuniformity then?
			int index = getColorIndex(pixel, true);
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
		maxbins = 0;

		for (int i = 0; i < 65536; ++i) {
			if (bins[i] == null)
				continue;

			double d = 1.0 / (double)bins[i].cnt;
			bins[i].ac *= d;
			bins[i].rc *= d;
			bins[i].gc *= d;
			bins[i].bc *= d;
			if (quan_sqrt)
				bins[i].cnt = (int) Math.sqrt(bins[i].cnt);
			bins[maxbins++] = bins[i];
		}

		for (int i = 0; i < maxbins - 1; i++) {
			bins[i].fw = (i + 1);
			bins[i + 1].bk = i;
		}
		// !!! Already zeroed out by calloc()
		//	bins[0].bk = bins[i].fw = 0;

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
		extbins = maxbins - nMaxColors;
		for (int i = 0; i < extbins; ) {
			/* Use heap to find which bins to merge */
			for (;;) {
				Pnnbin tb = bins[b1 = heap[1]]; /* One with least error */
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
			Pnnbin tb = bins[b1];
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
		List<Color> palette = new ArrayList<>();
		short k = 0;
		for (int i = 0;; ++k) {
			int alpha = (int) Math.rint(bins[i].ac);
			palette.add(new Color((int) Math.rint(bins[i].rc), (int) Math.rint(bins[i].gc), (int) Math.rint(bins[i].bc), alpha));
			if (hasTransparency && palette.get(k).equals(m_transparentColor)) {
				Color temp = palette.get(0);
				palette.set(0, palette.get(k));
				palette.set(k, temp);
			}

			if ((i = bins[i].fw) == 0)
				break;
		}

		return palette.toArray(new Color[0]);
	}
	
	private  int colorIndexToRGBA(final Color[] palette, final int k)
	{
		Color c1 = palette[k];
		return (c1.getAlpha() << 24) | (c1.getRed() << 16) | (c1.getGreen() << 8) | c1.getBlue();
	}

	private int nearestColorIndex(final Color[] palette, final int[] squares3, final Color c)
	{
		int k = 0;
		int curdist, mindist = SHORT_MAX;
		for (int i=0; i<palette.length; ++i) {
			Color c2 = palette[i];
			
			int adist = Math.abs(c2.getAlpha() - c.getAlpha());
			curdist = squares3[adist];
			if (curdist > mindist)
				continue;

			int rdist = Math.abs(c2.getRed() - c.getRed());
			curdist += squares3[rdist];
			if (curdist > mindist)
				continue;

			int gdist = Math.abs(c2.getGreen() - c.getGreen());
			curdist += squares3[gdist];
			if (curdist > mindist)
				continue;

			int bdist = Math.abs(c2.getBlue() - c.getBlue());
			curdist += squares3[bdist];
			if (curdist > mindist)
				continue;

			mindist = curdist;
			k = i;
		}
		return k;
	}

	private int closestColorIndex(final Color[] palette, final int[] squares3, final Color c)
	{
		int k = 0;
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
					closest[0] = (short) k;
					closest[2] = closest[4];
				}
				else if (closest[4] < closest[3]) {
					closest[1] = (short) k;
					closest[3] = closest[4];
				}
			}

			if (closest[3] == SHORT_MAX)
				closest[2] = 0;
		}
		else
			closest = got;

		if (closest[2] == 0 || (Math.random() % (closest[3] + closest[2])) <= closest[3])
			k = closest[0];
		else
			k = closest[1];

		closestMap.put(c, closest);
		return k;
	}

	boolean quantize_image(final Color[] pixels, final Color[] palette, int[] qPixels, final int width, final int height, final boolean dither)
	{
		int nMaxColors = palette.length;
		int[] sqr_tbl = new int[BYTE_MAX + BYTE_MAX + 1];

		for (int i = (-BYTE_MAX); i <= BYTE_MAX; i++)
			sqr_tbl[i + BYTE_MAX] = i * i;

		int[] squares3 = new int[sqr_tbl.length - BYTE_MAX];
		for (int i = 0; i < squares3.length; i++)
			squares3[i] = sqr_tbl[i + BYTE_MAX];

		int pixelIndex = 0;
		if (dither) {
			boolean odd_scanline = false;
			short[] row0, row1;
			int a_pix, r_pix, g_pix, b_pix, dir, k;
			final int DJ = 4;
			final int DITHER_MAX = 20;
			final int err_len = (width + 2) * DJ;
			short[] clamp = new short[DJ * 256];
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
					int offset = getColorIndex(c1, true);
					if (lookup[offset] == 0)
						lookup[offset] = nearestColorIndex(palette, squares3, c1) + 1;
					qPixels[pixelIndex] = lookup[offset] - 1;

					Color c2 = palette[qPixels[pixelIndex]];

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
				qPixels[i] = nearestColorIndex(palette, squares3, pixels[i]);
		}
		else {
			for (int i = 0; i < qPixels.length; i++)
				qPixels[i] = closestColorIndex(palette, squares3, pixels[i]);
		}

		return true;
	}

	boolean quantize_image(final Color[] pixels, int[] qPixels, int width, int height)
	{
		int[] squares3 = new int[BYTE_MAX + 1];
		for (int i = 0; i < squares3.length; i++)
			squares3[i] = i * i;

		short pixelIndex = 0;
		boolean odd_scanline = false;
		short[] row0, row1;
		int a_pix, r_pix, g_pix, b_pix, dir, k;
		final int DJ = 4;
		final int DITHER_MAX = 20;
		final int err_len = (width + 2) * DJ;
		short[] clamp = new short[DJ * 256];
		int[] limtb = new int[512];
		short[] erowerr = new short[err_len];
		short[] orowerr = new short[err_len];
		Color[] lookup = new Color[65536];

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

		for (int i = 0; i < height; i++) {
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
				int offset = getColorIndex(c1, false);
				if (lookup[offset] == null) {
					Color rgba1 = new Color((c1.getRed() & 0xF8), (c1.getGreen() & 0xFC), (c1.getBlue() & 0xF8), BYTE_MAX);
					if (hasSemiTransparency)
						rgba1 = new Color((c1.getRed() & 0xF0), (c1.getGreen() & 0xF0), (c1.getBlue() & 0xF0), (c1.getAlpha() & 0xF0));
					else if (hasTransparency)
						rgba1 = new Color((c1.getRed() & 0xF8), (c1.getGreen() & 0xF8), (c1.getBlue() & 0xF8), (c1.getAlpha() < BYTE_MAX) ? 0 : BYTE_MAX);
					lookup[offset] = rgba1;
				}
				qPixels[pixelIndex] = offset;

				Color c2 = lookup[offset];

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

	public int[] convert (int w, int h, int nMaxColors, boolean dither) {
		final Color[] cPixels = new Color[pixels.length];		
		for (int i =0; i<cPixels.length; ++i) {
			int pixel = pixels[i];
			int alfa = (pixel >> 24) & 0xff;
			int r   = (pixel >> 16) & 0xff;
			int g = (pixel >>  8) & 0xff;
			int b  = (pixel      ) & 0xff;
			cPixels[i] = new Color(r, g, b, alfa);
			if (alfa < BYTE_MAX) {
				hasSemiTransparency = true;
				if (alfa == 0) {
					hasTransparency = true;
					m_transparentColor = cPixels[i];
				}
			}			
		}
		
		if (nMaxColors > 256) {
			int[] qPixels = new int[cPixels.length];		
			quantize_image(cPixels, qPixels, w, h);
			return qPixels;
		}
		
		Pnnbin[] bins = new Pnnbin[65536];
		boolean quan_sqrt = nMaxColors > BYTE_MAX;
		Color[] palette = new Color[nMaxColors];
		if (nMaxColors > 2)
			palette = pnnquan(cPixels, bins, nMaxColors, quan_sqrt);
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

		int[] qPixels = new int[cPixels.length];		
		quantize_image(cPixels, palette, qPixels, w, h, true);
		closestMap.clear();
		
		for (int i = 0; i < qPixels.length; i++)
			qPixels[i] = colorIndexToRGBA(palette, qPixels[i]);
		
		return qPixels;
	}

}
