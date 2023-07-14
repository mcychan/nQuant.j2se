package nQuant.j2se;
/* Generalized Hilbert ("gilbert") space-filling curve for rectangular domains of arbitrary (non-power of two) sizes.
Copyright (c) 2021 - 2023 Miller Cy Chan
* A general rectangle with a known orientation is split into three regions ("up", "right", "down"), for which the function calls itself recursively, until a trivial path can be produced. */

import java.awt.Color;
import java.util.ArrayDeque;
import java.util.Queue;

public class GilbertCurve {

	private static final class ErrorBox
	{
		private final float[] p;
		private ErrorBox() {
			p = new float[4];
		}
		
		private ErrorBox(Color c) {
			p = new float[] {
				c.getRed(),
				c.getGreen(),
				c.getBlue(),
				c.getAlpha()
			};
		}
	}
	
	private byte ditherEdge;
	private final int width;
	private final int height;
	private final int[] pixels;
	private final Color[] palette;
	private final short[] qPixels;
	private final Ditherable ditherable;
	private final float[] saliencies;
	private final Queue<ErrorBox> errorq;
	private final float[] weights;
	private final int[] lookup;	
	
	private final byte DITHER_MAX, ditherMax;	
	private final int thresold;
	private static final float BLOCK_SIZE = 343f;
	
	private GilbertCurve(final int width, final int height, final int[] pixels, final Color[] palette, final short[] qPixels, final Ditherable ditherable, final float[] saliencies, double weight)
	{
		this.width = width;
		this.height = height;
		this.pixels = pixels;
		this.palette = palette;
		this.qPixels = qPixels;
		this.ditherable = ditherable;
		this.saliencies = saliencies;
		errorq = new ArrayDeque<>();
		boolean hasAlpha = weight < 0;
		weight = Math.abs(weight);
		DITHER_MAX = weight < .01 ? (weight > .0025) ? (byte) 25 : 16 : 9;
		double edge = hasAlpha ? 1 : Math.exp(weight) - .25;
		ditherMax = ditherEdge = (hasAlpha || DITHER_MAX > 9) ? (byte) BitmapUtilities.sqr(Math.sqrt(DITHER_MAX) + edge) : DITHER_MAX;
		if(palette.length / weight > 5000 && (weight > .045 || (weight > .01 && palette.length <= 64)))
			ditherEdge = (byte) BitmapUtilities.sqr(5 + edge);
		thresold = DITHER_MAX > 9 ? -112 : -88;
		weights = new float[DITHER_MAX];
		lookup = new int[65536];
	}

	private void ditherPixel(int x, int y){
		final int bidx = x + y * width;
		Color pixel = new Color(pixels[bidx], true);
		ErrorBox error = new ErrorBox(pixel);
		int i = 0;
		float maxErr = DITHER_MAX - 1; 
		for(ErrorBox eb : errorq) { 
			for(int j = 0; j < eb.p.length; ++j) {
				error.p[j] += eb.p[j] * weights[i];
				if(error.p[j] > maxErr)
					maxErr = error.p[j];
			}
			++i;
		}

		int r_pix = (int) Math.min(0xFF, Math.max(error.p[0], 0.0));
		int g_pix = (int) Math.min(0xFF, Math.max(error.p[1], 0.0));
		int b_pix = (int) Math.min(0xFF, Math.max(error.p[2], 0.0));
		int a_pix = (int) Math.min(0xFF, Math.max(error.p[3], 0.0));
		
		Color c2 = new Color(r_pix, g_pix, b_pix, a_pix);
		if (palette.length <= 32 && a_pix > 0xF0) {
			int offset = ditherable.getColorIndex(c2);
			if (lookup[offset] == 0)
				lookup[offset] = ditherable.nearestColorIndex(palette, c2, bidx) + 1;
			qPixels[bidx] = (short) (lookup[offset] - 1);
			
			if(saliencies != null && saliencies[bidx] > .65f && saliencies[bidx] < .75f) {
				final float strength = 1 / 3f;
				c2 = BlueNoise.diffuse(pixel, palette[qPixels[bidx]], 1 / saliencies[bidx], strength, x, y);
				qPixels[bidx] = ditherable.nearestColorIndex(palette, c2, bidx);
			}
		}
		else
			qPixels[bidx] = ditherable.nearestColorIndex(palette, c2, bidx);		
		
		errorq.poll();
		c2 = palette[qPixels[bidx]];
		if (palette.length > 256)
			qPixels[bidx] = (short) ditherable.getColorIndex(c2);
		
		error.p[0] = r_pix - c2.getRed();
		error.p[1] = g_pix - c2.getGreen();
		error.p[2] = b_pix - c2.getBlue();
		error.p[3] = a_pix - c2.getAlpha();

		boolean denoise = palette.length > 2;		
		boolean diffuse = BlueNoise.RAW_BLUE_NOISE[bidx & 4095] > thresold;
		double yDiff = diffuse ? 1 : CIELABConvertor.Y_Diff(pixel, c2);
		boolean illusion = !diffuse && BlueNoise.RAW_BLUE_NOISE[(int) (yDiff * 4096)] > thresold;
		byte ditherCursor = illusion ? ditherMax : ditherEdge;
		
		int errLength = denoise ? error.p.length - 1 : 0;	
		for(int j = 0; j < errLength; ++j) {
			if(Math.abs(error.p[j]) >= ditherCursor) {
				if (diffuse)
					error.p[j] = (float) Math.tanh(error.p[j] / maxErr * 20) * (ditherCursor - 1);
				else {					
					if(illusion)
						error.p[j] /= (float) (1 + Math.sqrt(ditherCursor));
					else
						error.p[j] = (float) (error.p[j] / maxErr * yDiff) * (ditherCursor - 1);						
				}
			}
		}

		errorq.add(error);
	}

	private void generate2d(int x, int y, int ax, int ay, int bx, int by) {
		int w = Math.abs(ax + ay);
		int h = Math.abs(bx + by);
		int dax = Integer.signum(ax);
		int day = Integer.signum(ay);
		int dbx = Integer.signum(bx);
		int dby = Integer.signum(by);

		if (h == 1) {
			for (int i = 0; i < w; ++i){
				ditherPixel(x, y);
				x += dax;
				y += day;
			}
			return;
		}

		if (w == 1) {
			for (int i = 0; i < h; ++i){
				ditherPixel(x, y);
				x += dbx;
				y += dby;
			}
			return;
		}

		int ax2 = ax / 2;
		int ay2 = ay / 2;
		int bx2 = bx / 2;
		int by2 = by / 2;

		int w2 = Math.abs(ax2 + ay2);
		int h2 = Math.abs(bx2 + by2);

		if (2 * w > 3 * h) {
			if ((w2 % 2) != 0 && w > 2) {
				ax2 += dax;
				ay2 += day;
			}
			generate2d(x, y, ax2, ay2, bx, by);
			generate2d(x + ax2, y + ay2, ax - ax2, ay - ay2, bx, by);
			return;
		}
		
		if ((h2 % 2) != 0 && h > 2) {
			bx2 += dbx;
			by2 += dby;
		}
		
		generate2d(x, y, bx2, by2, ax2, ay2);
		generate2d(x + bx2, y + by2, ax, ay, bx - bx2, by - by2);
		generate2d(x + (ax - dax) + (bx2 - dbx), y + (ay - day) + (by2 - dby), -bx2, -by2, -(ax - ax2), -(ay - ay2));
	}

	private void run() {
		/* Dithers all pixels of the image in sequence using
		 * the Gilbert path, and distributes the error in
		 * a sequence of DITHER_MAX pixels.
		 */
		final float weightRatio = (float) Math.pow(BLOCK_SIZE + 1f, 1f / (DITHER_MAX - 1f));
		float weight = 1f, sumweight = 0f;
		for(int c = 0; c < DITHER_MAX; ++c) {
			errorq.add(new ErrorBox());
			sumweight += (weights[DITHER_MAX - c - 1] = weight);
			weight /= weightRatio;
		}
		
		weight = 0f; /* Normalize */
		for(int c = 0; c < DITHER_MAX; ++c)
			weight += (weights[c] /= sumweight);
		weights[0] += 1f - weight;
		
		if (width >= height)
			generate2d(0, 0, width, 0, 0, height);
		else
			generate2d(0, 0, 0, height, width, 0);
	}
	
	public static short[] dither(final int width, final int height, final int[] pixels, final Color[] palette, final Ditherable ditherable, final float[] saliencies, final double weight)
	{
		short[] qPixels = new short[pixels.length];
		new GilbertCurve(width, height, pixels, palette, qPixels, ditherable, saliencies, weight).run();
		return qPixels;
	}
}
