package nQuant.j2se;

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.awt.image.DataBuffer;
import java.awt.image.IndexColorModel;
import java.util.HashMap;
import java.util.Map;

import static nQuant.j2se.BitmapUtilities.BYTE_MAX;

public class Otsu
{
	protected short alphaThreshold = 0;
	protected boolean hasSemiTransparency = false;
	protected int m_transparentPixelIndex = -1;
	protected Color m_transparentColor;
	
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
	private void getHistogram(Color[] pixels, int[] hist)
	{
		for(Color c : pixels)
		{
			if(c.getAlpha() <= alphaThreshold)
				continue;
			
			hist[c.getRed()]++;
			hist[c.getGreen()]++;
			hist[c.getBlue()]++;
		}
	}

	private short getOtsuThreshold(Color[] pixels)
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

	private boolean threshold(Color[] pixels, short thresh, float weight)
	{
		int maxThresh = (int) thresh;
		if (thresh >= 200)
		{
			weight = .85f;
			maxThresh = (int) (thresh * weight);
			thresh = 200;				
		}

		int minThresh = (int)(thresh * weight);		
		for (int i = 0; i < pixels.length; ++i)
		{
			Color c = pixels[i];
			if (c.getRed() + c.getGreen() + c.getBlue() > maxThresh * 3)
				pixels[i] = new Color(BYTE_MAX, BYTE_MAX, BYTE_MAX, c.getAlpha());
			else if (m_transparentPixelIndex >= 0 || c.getRed() + c.getGreen() + c.getBlue() < minThresh * 3)
				pixels[i] = new Color(0, 0, 0, c.getAlpha());
		}

		return true;
	}
	
	private boolean threshold(Color[] pixels, short thresh)
	{
		return threshold(pixels, thresh, 1.0f);
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
	
	private Color[] grabPixels(int[] pixels)
	{
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
		return cPixels;
	}

	public BufferedImage convertToGrayScale(BufferedImage srcimg)
	{
		int iWidth = srcimg.getWidth();
		int iHeight = srcimg.getHeight();
		int[] pixels = srcimg.getRGB(0, 0, iWidth, iHeight, null, 0, iWidth);

		float min1 = BYTE_MAX;
		float max1 = .0f;
		final Color[] cPixels = new Color[pixels.length];		
		for (int i = pixels.length - 1; i >= 0; --i) {
			int pixel = pixels[i];
			int alfa = (pixel >> 24) & 0xff;
			int green = (pixel >> 8) & 0xff;
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
			
			if(alfa <= alphaThreshold)
				continue;
			
			if (min1 > green)
				min1 = green;

			if (max1 < green)
				max1 = green;
		}

		for (int i = 0; i < pixels.length; ++i) {
			int green = cPixels[i].getGreen();
			int pixel = (int)((green - min1) * (BYTE_MAX / (max1 - min1)));
			pixels[i] = new Color(pixel, pixel, pixel, cPixels[i].getAlpha()).getRGB();
		}
		
		BufferedImage grayScaleImage = new BufferedImage(iWidth, iHeight, BufferedImage.TYPE_USHORT_GRAY);
		grayScaleImage.getRaster().setPixels(0, 0, iWidth, iHeight, pixels);
		return grayScaleImage;
	}	

	protected int getColorIndex(final Color c)
	{
		return BitmapUtilities.getColorIndex(c, hasSemiTransparency, m_transparentPixelIndex > -1);
	}
	
	protected Ditherable getDitherFn() {
		return new Ditherable() {
			@Override
			public int getColorIndex(Color c) {
				return BitmapUtilities.getColorIndex(c, hasSemiTransparency, m_transparentPixelIndex >= 0);
			}
			
			@Override
			public short nearestColorIndex(Color[] palette, Color c) {
				return Otsu.this.nearestColorIndex(palette, c);
			}
			
		};
	}

	public BufferedImage convertGrayScaleToBinary(BufferedImage srcimg, boolean isGrayscale)
	{
		BufferedImage sourceImg = isGrayscale ? srcimg : convertToGrayScale(srcimg);						

		int bitmapWidth = sourceImg.getWidth();
		int bitmapHeight = sourceImg.getHeight();

		int[] pixels = srcimg.getRGB(0, 0, bitmapWidth, bitmapHeight, null, 0, bitmapWidth);
		Color[] cPixels = grabPixels(pixels);

		short otsuThreshold = getOtsuThreshold(cPixels);
		if (!threshold(cPixels, otsuThreshold))
			return sourceImg;

		Color[] palette = new Color[2];
		if (m_transparentPixelIndex >= 0)
		{
			palette[0] = m_transparentColor;
			palette[1] = Color.BLACK;
		}
		else {
			palette[0] = Color.BLACK;
			palette[1] = Color.WHITE;
		}
		setColorModel(palette);	

		short[] qPixels = GilbertCurve.dither(bitmapWidth, bitmapHeight, cPixels, palette, getDitherFn());
		if (m_transparentPixelIndex >= 0) {
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
