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

	private void threshold(int[] pixels, short thresh, float weight)
	{
		int maxThresh = (int) thresh;
		if (thresh >= 200)
		{
			weight = hasAlpha() ? .9f : .8f;
			maxThresh = (int) (thresh * weight);
			thresh = 200;
		}

		int minThresh = (int)(thresh * weight);		
		for (int i = 0; i < pixels.length; ++i)
		{
			Color c = new Color(pixels[i], true);
			if (c.getRed() + c.getGreen() + c.getBlue() > maxThresh * 3)
				pixels[i] = new Color(BYTE_MAX, BYTE_MAX, BYTE_MAX, c.getAlpha()).getRGB();
			else if (hasAlpha() || c.getRed() + c.getGreen() + c.getBlue() < minThresh * 3)
				pixels[i] = new Color(0, 0, 0, c.getAlpha()).getRGB();
		}
	}
	
	private void threshold(int[] pixels, short thresh)
	{
		threshold(pixels, thresh, 1.0f);
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
	
	private void convertToGrayScale(int[] pixels)
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
			pixels[i] = new Color(grey, grey, grey, alfa).getRGB();
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
		
		if(!isGrayscale)
			convertToGrayScale(pixels);

		short otsuThreshold = getOtsuThreshold(pixels);
		threshold(pixels, otsuThreshold);

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
