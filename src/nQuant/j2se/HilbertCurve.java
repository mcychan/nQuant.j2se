package nQuant.j2se;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;

import static nQuant.j2se.HilbertCurve.Direction.*;

public class HilbertCurve {
	public enum Direction { LEFT, RIGHT, DOWN, UP };
	
	private final class ErrorBox
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
	
	private int x, y;
	private final int width;
	private final int height;
	private final Color[] image;
	private final Color[] palette;
	private final short[] qPixels;
	private final Ditherable ditherable;
	private final List<ErrorBox> errorq;
	private final float[] weights;
    
	private static final byte DITHER_MAX = 16;
	private static final float BLOCK_SIZE = 256f;	    
    
    private HilbertCurve(final int width, final int height, final Color[] image, final Color[] palette, final short[] qPixels, final Ditherable ditherable)
    {
    	x = 0;
    	y = 0;
    	this.width = width;
    	this.height = height;
        this.image = image;
        this.palette = palette;
        this.qPixels = qPixels;
        this.ditherable = ditherable;	        
        errorq = new ArrayList<>();
        weights = new float[DITHER_MAX];
    }
    
    private void ditherCurrentPixel()
	{
	    if(x >= 0 && y >= 0 && x < width && y < height) {
	    	Color pixel = image[x + y * width];
	    	ErrorBox error = new ErrorBox(pixel);	    	
	        for(int c = 0; c < DITHER_MAX; ++c) {
	        	ErrorBox eb = errorq.get(c);
	        	for(int j = 0; j < eb.p.length; ++j)
	        		error.p[j] += eb.p[j] * weights[c];
	        }

	        int r_pix = (int) Math.min(BLOCK_SIZE-1, Math.max(error.p[0], 0.0));
	        int g_pix = (int) Math.min(BLOCK_SIZE-1, Math.max(error.p[1], 0.0));
	        int b_pix = (int) Math.min(BLOCK_SIZE-1, Math.max(error.p[2], 0.0));
	        int a_pix = (int) Math.min(BLOCK_SIZE-1, Math.max(error.p[3], 0.0));
	        
	        Color c2 = new Color(r_pix, g_pix, b_pix, a_pix);		        
	        qPixels[x + y * width] = ditherable.nearestColorIndex(palette, c2);

	        errorq.remove(0);
	        c2 = palette[qPixels[x + y * width]];
	        error.p[0] = r_pix - c2.getRed();
	        if(Math.abs(error.p[0]) > DITHER_MAX)
	        	error.p[0] *= .8f;
	        error.p[1] = g_pix - c2.getGreen();
	        if(Math.abs(error.p[1]) > DITHER_MAX)
	        	error.p[1] *= .8f;
	        error.p[2] = b_pix - c2.getBlue();
	        if(Math.abs(error.p[2]) > DITHER_MAX)
	        	error.p[2] *= .8f;
	        error.p[3] = a_pix - c2.getAlpha();
	        if(Math.abs(error.p[3]) > DITHER_MAX)
	        	error.p[3] *= .8f;
	        errorq.add(error);
	    }
	}
    
    private void run()
    {
        /* Dithers all pixels of the image in sequence using
         * the Hilbert path, and distributes the error in
         * a sequence of 16 pixels.
         */
        x = y = 0;
        final float weightRatio = (float) Math.pow(BLOCK_SIZE + 1f,  1f / (DITHER_MAX - 1f));
        float weight = 1f, sumweight = 0f;
        for(int c = 0; c < DITHER_MAX; ++c)
        {
            errorq.add(new ErrorBox());
            sumweight += (weights[DITHER_MAX - c - 1] = 1.0f / weight);
            weight *= weightRatio;
        }
        
        weight = 0f; /* Normalize */
        for(int c = 0; c < DITHER_MAX; ++c)
            weight += (weights[c] /= sumweight);
        weights[0] += 1f - weight;
        /* Walk the path. */
        int i = Math.max(width, height), depth = 0;
        while(i > 0) {
        	++depth;
        	i >>= 1;
        }
        
        iter(depth, UP);
        ditherCurrentPixel();
    }
    
    private void navTo(Direction dir)
	{
    	ditherCurrentPixel();
		switch(dir)
        {
            case LEFT:
            	--x;
            	break;
            case RIGHT:
            	++x;
            	break;
            case UP:
            	--y;
            	break;
            case DOWN:
            	++y;
            	break;
        }
	}
	
    private void curve(final int level, Direction a, Direction b, Direction c, Direction d, Direction e, Direction f, Direction g)
    {
		iter(level-1, a);
		navTo(e);
		iter(level-1, b);
		navTo(f);
        iter(level-1, c);
        navTo(g);
        iter(level-1, d);
    }

    private void iter(final int level, Direction dir) {
    	if(level <= 0)
    		return;    	

        switch(dir)
        {
            case LEFT:
            	curve(level, UP,LEFT,LEFT,DOWN, RIGHT,DOWN,LEFT);
            	break;
            case RIGHT:
            	curve(level, DOWN,RIGHT,RIGHT,UP, LEFT,UP,RIGHT);
            	break;
            case UP:
            	curve(level, LEFT,UP,UP,RIGHT, DOWN,RIGHT,UP);
            	break;
            case DOWN:
            	curve(level, RIGHT,DOWN,DOWN,LEFT, UP,LEFT,DOWN);
            	break;
        }
    }
    
    public static short[] dither(final int width, final int height, final Color[] pixels, final Color[] palette, final Ditherable ditherable)
    {
    	short[] qPixels = new short[pixels.length];		
    	new HilbertCurve(width, height, pixels, palette, qPixels, ditherable).run();        
        return qPixels;
    }
}
