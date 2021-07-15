package nQuant.j2se;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;

public class HilbertCurve {
	private enum Direction { W,E,S,N };
	
	private class ErrorBox
	{
		double p[] = new double[4];
	}
	
	private int x, y;
	private final int width;
	private final int height;
	private final Color[] image;
	private final Color[] palette;
	private final short[] qPixels;
	private final Ditherable ditherable;
	private final List<ErrorBox> errorq;
	private final double[] weights;
    
	private static final byte DITHER_MAX = 16;
	private static final double BLOCK_SIZE = 256.0;	    
    
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
        weights = new double[DITHER_MAX];
    }	    
    
    private void dither_thispix()
	{
	    if(x >= 0 && y >= 0 && x < width && y < height) {
	    	ErrorBox error = new ErrorBox();
	    	Color pixel = image[x + y * width];
	    	error.p[0] = pixel.getRed();
	    	error.p[1] = pixel.getGreen();
	    	error.p[2] = pixel.getBlue();
	    	error.p[3] = pixel.getAlpha();
	    	
	        for(int c = 0; c < DITHER_MAX; ++c) {
	        	for(int j = 0; j < 4; ++j)
	        		error.p[j] += errorq.get(c).p[j] * weights[c];
	        }

	        int r = (int)Math.min(BLOCK_SIZE-1, Math.max(error.p[0], 0.0));
	        int g = (int)Math.min(BLOCK_SIZE-1, Math.max(error.p[1], 0.0));
	        int b = (int)Math.min(BLOCK_SIZE-1, Math.max(error.p[2], 0.0));
	        int a = (int)Math.min(BLOCK_SIZE-1, Math.max(error.p[3], 0.0));
	        Color c2 = new Color(r, g, b, a);		        
	        qPixels[x + y * width] = ditherable.nearestColorIndex(palette, c2);

	        errorq.remove(0);
	        c2 = palette[qPixels[x + y * width]];
	        error.p[0] = r - c2.getRed();
	        error.p[1] = g - c2.getGreen();
	        error.p[2] = b - c2.getBlue();
	        error.p[3] = a - c2.getAlpha();
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
        final double weightRatio = Math.pow(BLOCK_SIZE + 1,  1 / (DITHER_MAX - 1.0));
        double weight = 1.0, sumweight = 0.0;
        for(int c = 0; c < DITHER_MAX; ++c)
        {
            errorq.add(new ErrorBox());
            sumweight += (weights[DITHER_MAX - c - 1] = 1.0 / weight);
            weight *= weightRatio;
        }
        weight = 0; /* Normalize */
        for(int c = 0; c < DITHER_MAX; ++c)
            weight += (weights[c] /= sumweight);
        weights[0] += 1.0 - weight;
        /* Walk the path. */
        int i = Math.max(width, height), depth = 0;
        while(i != 0) {
        	++depth;
        	i /= 2;
        }
        
        iter(depth, Direction.N);
        dither_thispix();
    }
    
    private void navTo(Direction dir)
	{
    	dither_thispix();
		switch(dir)
        {
            case W: --x;
            	break;
            case E: ++x;
            	break;
            case N: --y;
            	break;
            case S: ++y;
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
            case W: curve(level, Direction.N,Direction.W,Direction.W,Direction.S, Direction.E,Direction.S,Direction.W);
            	break;
            case E: curve(level, Direction.S,Direction.E,Direction.E,Direction.N, Direction.W,Direction.N,Direction.E);
            	break;
            case N: curve(level, Direction.W,Direction.N,Direction.N,Direction.E, Direction.S,Direction.E,Direction.N);
            	break;
            case S: curve(level, Direction.E,Direction.S,Direction.S,Direction.W, Direction.N,Direction.W,Direction.S);
            	break;
        }
    }
    
    public static short[] dither(final int width, final int height, final Color[] pixels, final Color[] palette, final Ditherable ditherable)
    {
    	short[] qPixels = new short[pixels.length];		
    	HilbertCurve hc = new HilbertCurve(width, height, pixels, palette, qPixels, ditherable);
    	hc.run();        
        return qPixels;
    }
}
