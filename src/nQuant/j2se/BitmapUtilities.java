package nQuant.j2se;

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.awt.image.ColorModel;
import java.awt.image.DataBufferByte;
import java.awt.image.DataBufferShort;
import java.awt.image.DirectColorModel;
import java.awt.image.IndexColorModel;
import java.awt.image.Raster;
import java.awt.image.WritableRaster;

public class BitmapUtilities {
	protected static final char BYTE_MAX = -Byte.MIN_VALUE + Byte.MAX_VALUE;

	public static double sqr(double value)
	{
		return value * value;
	}
	
	public static int getBitsPerPixel(int nMaxColors)
	{
		return (int) Math.ceil(Math.log(nMaxColors) / Math.log(2));
	}
	
	protected static int[] calcDitherPixel(Color c, short[] clamp, int[] rowerr, int cursor, boolean noBias)
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
	
	protected static short[] quantize_image(final int width, final int height, final int[] pixels, final Color[] palette, final Ditherable ditherable, final boolean hasSemiTransparency, final boolean dither)
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
					Color c = new Color(pixels[pixelIndex], true);
					int[] ditherPixel = calcDitherPixel(c, clamp, row0, cursor0, noBias);
					int r_pix = ditherPixel[0];
				int g_pix = ditherPixel[1];
				int b_pix = ditherPixel[2];
				int a_pix = ditherPixel[3];

					Color c1 = new Color(r_pix, g_pix, b_pix, a_pix);
					if(noBias && a_pix > 0xF0) {
						int offset = ditherable.getColorIndex(c1);
						if (lookup[offset] == 0)
							lookup[offset] = (c.getAlpha() == 0) ? 1 : ditherable.nearestColorIndex(palette, c1, i + j) + 1;
						qPixels[pixelIndex] = (short) (lookup[offset] - 1);
					}
					else
						qPixels[pixelIndex] = (c.getAlpha() == 0) ? 0 : ditherable.nearestColorIndex(palette, c1, i + j);

					Color c2 = palette[qPixels[pixelIndex]];
					if(nMaxColors > 256)
						qPixels[pixelIndex] = (short) ditherable.getColorIndex(c2);

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

		for (int i = 0; i < qPixels.length; ++i)
			qPixels[i] = ditherable.nearestColorIndex(palette, new Color(pixels[i], true), i);

		return qPixels;
	}
	
	public static int getColorIndex(final Color c, boolean hasSemiTransparency, boolean hasTransparency)
	{
		if(hasSemiTransparency)
			return (c.getAlpha() & 0xF0) << 8 | (c.getRed() & 0xF0) << 4 | (c.getGreen() & 0xF0) | (c.getBlue() >> 4);
		if (hasTransparency)
			return (c.getAlpha() & 0x80) << 8 | (c.getRed() & 0xF8) << 7 | (c.getGreen() & 0xF8) << 2 | (c.getBlue() >> 3);
		return (c.getRed() & 0xF8) << 8 | (c.getGreen() & 0xFC) << 3 | (c.getBlue() >> 3);
	}
	
	public static BufferedImage toIndexedBufferedImage(short[] qPixels, IndexColorModel icm, int width, int height) {		
		if(icm.getPixelSize() < 8) {			
			WritableRaster raster = Raster.createWritableRaster(icm.createCompatibleSampleModel(width, height), new DataBufferByte(qPixels.length), null);
			for(int y = 0; y < height; ++y) {
		    	for(int x = 0; x < width; ++x)
		    		raster.setSample(x, y, 0, qPixels[x + y * width]);
			}
			return new BufferedImage(icm, raster, icm.isAlphaPremultiplied(), null);
		}
		
		WritableRaster raster = Raster.createWritableRaster(icm.createCompatibleSampleModel(width, height), new DataBufferShort(qPixels, qPixels.length), null);
		return new BufferedImage(icm, raster, icm.isAlphaPremultiplied(), null);
	}
	
	public static BufferedImage processImagePixels(short[] qPixels, ColorModel colorModel, int width, int height) {
		if(colorModel instanceof IndexColorModel)
			return toIndexedBufferedImage(qPixels, (IndexColorModel) colorModel, width, height);

		BufferedImage highColorImage = null;
		if(colorModel instanceof DirectColorModel) {
			ColorModel cmSw = colorModel;
			WritableRaster wr = cmSw.createCompatibleWritableRaster(width, height);
			highColorImage = new BufferedImage(cmSw, wr, cmSw.isAlphaPremultiplied(), null);
		}
		else
			highColorImage = new BufferedImage(width, height, BufferedImage.TYPE_USHORT_565_RGB);

		highColorImage.getRaster().setDataElements(0, 0, width, height, qPixels);
		return highColorImage;
	}
}
