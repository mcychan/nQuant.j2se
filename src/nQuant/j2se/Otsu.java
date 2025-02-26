package nQuant.j2se;

/* Otsu's Image Segmentation Method
* Copyright (C) 2009 Tolga Birdal
* Copyright (c) 2018 - 2025 Miller Cy Chan
*/

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.awt.image.DataBuffer;
import java.awt.image.IndexColorModel;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import static nQuant.j2se.BitmapUtilities.BYTE_MAX;

public class Otsu
{
	protected short alphaThreshold = 0xF;
	protected boolean hasSemiTransparency = false;
	protected int m_transparentPixelIndex = -1;
	protected Color m_transparentColor = new Color(BYTE_MAX, BYTE_MAX, BYTE_MAX, 0);
	
	protected IndexColorModel m_colorModel;
	protected Map<Integer, Short> nearestMap = new HashMap<>();	
	
	protected final void setColorModel(final Color[] palette)
	{
		int nMaxColors = palette.length;
		
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

	// function is used to compute the q values in the equation
	private static float px(int init, int end, int[] hist)
	{
		int sum = 0;
		for (int i = init; i <= end; ++i)
			sum += hist[i];

		return sum;
	}

	// function is used to compute the mean values in the equation (mu)
	private static float mx(int init, int end, int[] hist)
	{
		int sum = 0;
		for (int i = init; i <= end; ++i)
			sum += i * hist[i];

		return sum;
	}

	// finds the maximum element in a vector
	private static short findMax(float[] vec, int n)
	{
		float maxVec = 0;
		short idx = 0;

		for (short i = 1; i < n - 1; ++i)
		{
			if (vec[i] > maxVec)
			{
				maxVec = vec[i];
				idx = i;
			}
		}
		return idx;
	}

	// simply computes the image histogram
	private void getHistogram(int[] pixels, int[] hist)
	{
		for(int pixel : pixels)
		{
			Color c = new Color(pixel, true);
			if(c.getAlpha() <= alphaThreshold)
				continue;
			
			hist[c.getRed()]++;
			hist[c.getGreen()]++;
			hist[c.getBlue()]++;
		}
	}

	private short getOtsuThreshold(int[] pixels)
	{
		float[] vet = new float[256];
		int[] hist = new int[256];

		getHistogram(pixels, hist);

		// loop through all possible t values and maximize between class variance
		for (int k = 1; k != BYTE_MAX; ++k)
		{
			float p1 = px(0, k, hist);
			float p2 = px(k + 1, BYTE_MAX, hist);
			float p12 = p1 * p2;
			if (p12 == 0)
				p12 = 1;
			float diff = (mx(0, k, hist) * p2) - (mx(k + 1, BYTE_MAX, hist) * p1);
			vet[k] = diff * diff / p12;
		}

		return findMax(vet, 256);
	}	

	private void threshold(final int[] pixels, int[] dest, short thresh, float weight)
	{
		int maxThresh = (int) thresh;
		if (thresh >= 200)
		{
			weight = .78f;
			maxThresh = (int) (thresh * weight);
			thresh = 200;
		}

		int minThresh = (int)(thresh * (m_transparentPixelIndex >= 0 ? .9f : weight));
		final double shadow = m_transparentPixelIndex >= 0 ? 3.5 : 3;
		for (int i = 0; i < pixels.length; ++i)
		{
			Color c = new Color(pixels[i], true);
			if (c.getAlpha() < alphaThreshold && c.getRed() + c.getGreen() + c.getBlue() > maxThresh * 3)
				dest[i] = new Color(BYTE_MAX, BYTE_MAX, BYTE_MAX, c.getAlpha()).getRGB();
			else if (c.getRed() + c.getGreen() + c.getBlue() < minThresh * shadow)
				dest[i] = new Color(0, 0, 0, c.getAlpha()).getRGB();
		}
	}
	
	private void threshold(final int[] pixels, int[] dest, short thresh)
	{
		threshold(pixels, dest, thresh, 1.0f);
	}
	
	private int[] cannyFilter(final int width, final int[] pixelsGray, double lowerThreshold, double higherThreshold) {
		final int height = pixelsGray.length / width;
		final int area = width * height;

		int[] pixelsCanny = new int[area];
		Arrays.fill(pixelsCanny, Color.WHITE.getRGB());

		int[][] gx = new int[][] {{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}};
		int[][] gy = new int[][] {{-1, -2, -1}, {0, 0, 0}, {1, 2, 1}};
		double[] G = new double[area];
		int[] theta = new int[area];
		double largestG = 0.0;

		// perform canny edge detection on everything but the edges
		for (int i = 1; i < height - 1; ++i) {
			for (int j = 1; j < width - 1; ++j) {
				// find gx and gy for each pixel
				double gxValue = 0.0;
				double gyValue = 0.0;
				for (int x = -1; x <= 1; ++x) {
					for (int y = -1; y <= 1; ++y) {
						final Color c = new Color(pixelsGray[(i + x) * width + j + y], true);
						gxValue += gx[1 - x][1 - y] * c.getGreen();
						gyValue += gy[1 - x][1 - y] * c.getGreen();
					}
				}

				final int center = i * width + j;
				// calculate G and theta
				G[center] = Math.sqrt(BitmapUtilities.sqr(gxValue) + BitmapUtilities.sqr(gyValue));
				double atanResult = Math.atan2(gyValue, gxValue) * 180.0 / Math.PI;
				theta[center] = (int)(180.0 + atanResult);

				if (G[center] > largestG)
					largestG = G[center];

				// setting the edges
				if (i == 1) {
					G[center - 1] = G[center];
					theta[center - 1] = theta[center];
				}
				else if (j == 1) {
					G[center - width] = G[center];
					theta[center - width] = theta[center];
				}
				else if (i == height - 1) {
					G[center + 1] = G[center];
					theta[center + 1] = theta[center];
				}
				else if (j == width - 1) {
					G[center + width] = G[center];
					theta[center + width] = theta[center];
				}

				// setting the corners
				if (i == 1 && j == 1) {
					G[center - width - 1] = G[center];
					theta[center - width - 1] = theta[center];
				}
				else if (i == 1 && j == width - 1) {
					G[center - width + 1] = G[center];
					theta[center - width + 1] = theta[center];
				}
				else if (i == height - 1 && j == 1) {
					G[center + width - 1] = G[center];
					theta[center + width - 1] = theta[center];
				}
				else if (i == height - 1 && j == width - 1) {
					G[center + width + 1] = G[center];
					theta[center + width + 1] = theta[center];
				}

				// to the nearest 45 degrees
				theta[center] = (int) Math.rint(theta[center] / 45) * 45;
			}
		}

		largestG *= .5;

		// non-maximum suppression
		for (int i = 1; i < height - 1; ++i) {
			for (int j = 1; j < width - 1; ++j) {
				final int center = i * width + j;
				if (theta[center] == 0 || theta[center] == 180) {
					if (G[center] < G[center - 1] || G[center] < G[center + 1])
						G[center] = 0;
				}
				else if (theta[center] == 45 || theta[center] == 225) {
					if (G[center] < G[center + width + 1] || G[center] < G[center - width - 1])
						G[center] = 0;
				}
				else if (theta[center] == 90 || theta[center] == 270) {
					if (G[center] < G[center + width] || G[center] < G[center - width])
						G[center] = 0;
				}
				else {
					if (G[center] < G[center + width - 1] || G[center] < G[center - width + 1])
						G[center] = 0;
				}

				int grey = BYTE_MAX - ((int) Math.rint(G[center] * (255.0 / largestG)) & BYTE_MAX);
				Color c = new Color(pixelsGray[center], true);
				pixelsCanny[center] = new Color(grey, grey, grey, c.getAlpha()).getRGB();
			}
		}

		int k = 0;
		double minThreshold = lowerThreshold * largestG, maxThreshold = higherThreshold * largestG;
		do {
			for (int i = 1; i < height - 1; ++i) {
				for (int j = 1; j < width - 1; ++j) {
					final int center = i * width + j;
					if (G[center] < minThreshold)
						G[center] = 0;
					else if (G[center] >= maxThreshold)
						continue;
					else if (G[center] < maxThreshold) {
						G[center] = 0;
						for (int x = -1; x <= 1; ++x) {
							for (int y = -1; y <= 1; y++) {
								if (x == 0 && y == 0)
									continue;
								if (G[center + x * width + y] >= maxThreshold) {
									G[center] = higherThreshold * largestG;
									k = 0;
									x = 2;
									break;
								}
							}
						}
					}
					
					int grey = BYTE_MAX - ((int) Math.rint(G[center] * (255.0 / largestG)) & BYTE_MAX);
					Color c = new Color(pixelsGray[center], true);
					pixelsCanny[center] = new Color(grey, grey, grey, c.getAlpha()).getRGB();
				}
			}
		} while (k++ < 100);
		return pixelsCanny;
	}

	protected short nearestColorIndex(Color[] palette, final Color c)
	{
		Short got = nearestMap.get(c.getRGB());
		if (got != null)
			return got;

		short k = 0;
		if (c.getAlpha() <= alphaThreshold)
            return k;

		double mindist = 1e100;
		for (int i = 0; i < palette.length; ++i)
		{
			Color c2 = palette[i];
			double curdist = BitmapUtilities.sqr(c2.getAlpha() - c.getAlpha());
			if (curdist > mindist)
				continue;

			curdist += BitmapUtilities.sqr(c2.getRed() - c.getRed());
			if (curdist > mindist)
				continue;

			curdist += BitmapUtilities.sqr(c2.getGreen() - c.getGreen());
			if (curdist > mindist)
				continue;

			curdist += BitmapUtilities.sqr(c2.getBlue() - c.getBlue());
			if (curdist > mindist)
				continue;

			mindist = curdist;
			k = (short) i;
		}
		nearestMap.put(c.getRGB(), k);
		return k;
	}
	
	private void grabPixels(int[] pixels)
	{	
		for (int i = pixels.length - 1; i >= 0; --i) {
			int pixel = pixels[i];
			int alfa = (pixel >> 24) & 0xff;
			pixels[i] = new Color(pixel, true).getRGB();
			if (alfa < 0xE0) {
				if (alfa == 0) {
					m_transparentPixelIndex = i;
					pixels[i] = m_transparentColor.getRGB();
				}
				else if(alfa > alphaThreshold)
					hasSemiTransparency = true;
			}			
		}
	}

	public BufferedImage convertToGrayScale(BufferedImage srcimg)
	{
		int iWidth = srcimg.getWidth();
		int iHeight = srcimg.getHeight();
		int[] pixels = srcimg.getRGB(0, 0, iWidth, iHeight, null, 0, iWidth);

		float min1 = BYTE_MAX;
		float max1 = .0f;
		for (int i = pixels.length - 1; i >= 0; --i) {
			int pixel = pixels[i];
			int alfa = (pixel >> 24) & 0xff;
			int green = (pixel >> 8) & 0xff;
			if (alfa < 0xE0) {
				if (alfa == 0) {
					m_transparentPixelIndex = i;
					m_transparentColor = new Color(pixels[i], true);
				}
				else if(alfa > 0xF)
					hasSemiTransparency = true;
			}		
			
			if(alfa <= alphaThreshold)
				continue;
			
			if (min1 > green)
				min1 = green;

			if (max1 < green)
				max1 = green;
		}

		for (int i = 0; i < pixels.length; ++i) {
			int alfa = (pixels[i] >> 24) & 0xff;
			if(alfa <= alphaThreshold)
				continue;
			
			int green = (pixels[i] >> 8) & 0xff;
			int pixel = (int)((green - min1) * (BYTE_MAX / (max1 - min1)));
			pixels[i] = new Color(pixel, pixel, pixel, alfa).getRGB();
		}
		
		BufferedImage grayScaleImage = new BufferedImage(iWidth, iHeight, BufferedImage.TYPE_USHORT_GRAY);
		grayScaleImage.getRaster().setPixels(0, 0, iWidth, iHeight, pixels);
		return grayScaleImage;
	}	
	
	private void convertToGrayScale(final int[] pixels, int[] dest)
	{
		float min1 = BYTE_MAX;
		float max1 = .0f;

		for (int pixel : pixels)
		{
			int alfa = (pixel >> 24) & 0xff;
			if (alfa <= alphaThreshold)
				continue;
			
			int green = (pixel >> 8) & 0xff;

			if (min1 > green)
				min1 = green;

			if (max1 < green)
				max1 = green;
		}

		for (int i = 0; i < pixels.length; ++i)
		{
			int alfa = (pixels[i] >> 24) & 0xff;
			if (alfa <= alphaThreshold)
				continue;
			
			int green = (pixels[i] >> 8) & 0xff;
			int grey = (int)((green - min1) * (BYTE_MAX / (max1 - min1)));
			dest[i] = new Color(grey, grey, grey, alfa).getRGB();
		}
	}

	protected int getColorIndex(final Color c)
	{
		return BitmapUtilities.getColorIndex(c, hasSemiTransparency, m_transparentPixelIndex > -1);
	}
	
	protected Ditherable getDitherFn() {
		return new Ditherable() {
			@Override
			public int getColorIndex(Color c) {
				return BitmapUtilities.getColorIndex(c, hasSemiTransparency, hasAlpha());
			}
			
			@Override
			public short nearestColorIndex(Color[] palette, Color c, final int pos) {
				return Otsu.this.nearestColorIndex(palette, c);
			}			
		};
	}

	public BufferedImage convertGrayScaleToBinary(BufferedImage srcimg, boolean isGrayscale)
	{	
		int bitmapWidth = srcimg.getWidth();
		int bitmapHeight = srcimg.getHeight();

		int[] pixels = srcimg.getRGB(0, 0, bitmapWidth, bitmapHeight, null, 0, bitmapWidth);
		grabPixels(pixels);
		
		int[] pixelsGray = pixels.clone();
		if(!isGrayscale)
			convertToGrayScale(pixels, pixelsGray);

		short otsuThreshold = getOtsuThreshold(pixelsGray);
		double lowerThreshold = 0.03, higherThreshold = 0.1;
		pixels = cannyFilter(bitmapWidth, pixelsGray, lowerThreshold, higherThreshold);
		threshold(pixelsGray, pixels, otsuThreshold);

		Color[] palette = new Color[2];
		if (hasAlpha())
		{
			palette[0] = m_transparentColor;
			palette[1] = Color.BLACK;
		}
		else {
			palette[0] = Color.BLACK;
			palette[1] = Color.WHITE;
		}
		setColorModel(palette);	

		short[] qPixels = GilbertCurve.dither(bitmapWidth, bitmapHeight, pixels, palette, getDitherFn(), null, 1);
		if (hasAlpha()) {
			short k = qPixels[m_transparentPixelIndex];
			if (!palette[k].equals(m_transparentColor)) {
				Color c1 = palette[0]; palette[0] = palette[1]; palette[1] = c1;
			}
		}	

		nearestMap.clear();
		return BitmapUtilities.toIndexedBufferedImage(qPixels, m_colorModel, bitmapWidth, bitmapHeight);
	}
	
	public BufferedImage convertGrayScaleToBinary(BufferedImage srcimg)
	{
		return convertGrayScaleToBinary(srcimg, false);
	}
	
	public boolean hasAlpha() {
		return m_transparentPixelIndex > -1;
	}
}
