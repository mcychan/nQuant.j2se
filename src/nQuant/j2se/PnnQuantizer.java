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

import static nQuant.j2se.BitmapUtilities.BYTE_MAX;

public class PnnQuantizer {	
	protected short alphaThreshold = 0;
	protected boolean hasSemiTransparency = false;
	protected int m_transparentPixelIndex = -1;
	protected final int width, height;	
	protected int[] pixels;
	protected Color m_transparentColor;
	
	private double PR = .2126, PG = .7152, PB = .0722;
	private Color[] m_palette;
	private ColorModel m_colorModel;
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
		for (int i = bin1.fw; i != 0; i = bins[i].fw) {
			double nerr = PR * BitmapUtilities.sqr(bins[i].rc - wr) + PG * BitmapUtilities.sqr(bins[i].gc - wg) + PB * BitmapUtilities.sqr(bins[i].bc - wb);
			if(hasSemiTransparency)
				nerr += BitmapUtilities.sqr(bins[i].ac - wa);
			
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
			
			m_colorModel = new IndexColorModel(BitmapUtilities.getBitsPerPixel(nMaxColors),         // bits per pixel
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
			
			int index = BitmapUtilities.getColorIndex(c, hasSemiTransparency, nMaxColors < 64 || m_transparentPixelIndex >= 0);
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
		if (BitmapUtilities.sqr(nMaxColors) / maxbins < .03)
			quan_rt = 0;		

		int j = 0;
		for (; j < maxbins - 1; ++j) {
			bins[j].fw = j + 1;
			bins[j + 1].bk = j;
			
			if (quan_rt > 0) {
				if(nMaxColors < 64)
					bins[j].cnt = (float) Math.sqrt(bins[j].cnt);
				else
					bins[j].cnt = (int) Math.sqrt(bins[j].cnt);
			}
			else if (quan_rt < 0)
				bins[j].cnt = (int) Math.cbrt(bins[j].cnt);
		}
		if (quan_rt > 0) {
			if(nMaxColors < 64)
				bins[j].cnt = (float) Math.sqrt(bins[j].cnt);
			else
				bins[j].cnt = (int) Math.sqrt(bins[j].cnt);
		}		
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
			tb.ac = d * Math.round(n1 * tb.ac + n2 * nb.ac);
			tb.rc = d * Math.round(n1 * tb.rc + n2 * nb.rc);
			tb.gc = d * Math.round(n1 * tb.gc + n2 * nb.gc);
			tb.bc = d * Math.round(n1 * tb.bc + n2 * nb.bc);
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
			if (m_transparentPixelIndex >= 0 && alpha == 0) {
				Color temp = palette[0]; palette[0] = m_transparentColor; palette[k] = temp;
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
		int[] closest = new int[4];
		int[] got = closestMap.get(c.getRGB());
		if (got == null) {
			closest[2] = closest[3] = Integer.MAX_VALUE;

			for (; k < palette.length; ++k) {
				Color c2 = palette[k];
				if(c2 == null)
					break;

				final int ERR = (int) (Math.abs(c.getAlpha() - c2.getAlpha()) + Math.abs(c.getRed() - c2.getRed()) + Math.abs(c.getGreen() - c2.getGreen()) + Math.abs(c.getBlue() - c2.getBlue()));
				if (ERR < closest[2]) {
					closest[1] = closest[0];
					closest[3] = closest[2];
					closest[0] = k;
					closest[2] = ERR;
				}
				else if (ERR < closest[3]) {
					closest[1] = k;
					closest[3] = ERR;
				}
			}

			if (closest[3] == Integer.MAX_VALUE)
				closest[1] = closest[0];
			
			closestMap.put(c.getRGB(), closest);
		}
		else
			closest = got;

		Random rand = new Random();
		if (closest[2] == 0 || (rand.nextInt(32767) % (closest[3] + closest[2])) <= closest[3]) {
			if(closest[2] > palette.length)
				return nearestColorIndex(palette, c);
			return (short) closest[0];
		}
		
		if(closest[3] > palette.length)
			return nearestColorIndex(palette, c);
		return (short) closest[1];
	}		
	
	protected Ditherable getDitherFn(final boolean dither) {
		return new Ditherable() {
			@Override
			public int getColorIndex(Color c) {
				return BitmapUtilities.getColorIndex(c, hasSemiTransparency, m_transparentPixelIndex >= 0);
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
		if(hasSemiTransparency)
			qPixels = GilbertCurve.dither(width, height, cPixels, palette, getDitherFn(dither), 1.25f);
		else if (nMaxColors < 64 && nMaxColors > 32)
			qPixels = BitmapUtilities.quantize_image(width, height, cPixels, palette, getDitherFn(dither), hasSemiTransparency, dither);
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

	public BufferedImage convert(int nMaxColors, boolean dither) {
		final Color[] cPixels = new Color[pixels.length];		
		for (int i = pixels.length - 1; i >= 0; --i) {
			int pixel = pixels[i];
			int alfa = (pixel >> 24) & 0xff;
			cPixels[i] = new Color(pixel, true);
			if (alfa < BYTE_MAX) {				
				if (alfa == 0) {
					m_transparentPixelIndex = i;
					m_transparentColor = cPixels[i];
					if(m_transparentColor.getRGB() < BYTE_MAX)
						cPixels[i] = m_transparentColor = new Color(51, 102, 102, alfa);
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

		return BitmapUtilities.processImagePixels(qPixels, m_colorModel, width, height);
	}

	public Color[] getPalette() {
		return m_palette;
	}
	
	public boolean hasAlpha() {
		return m_transparentPixelIndex > -1;
	}

}
