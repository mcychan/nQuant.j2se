package nQuant.j2se;

import java.awt.Canvas;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.RenderingHints;
import java.awt.Toolkit;
import java.awt.image.BufferedImage;
import java.awt.image.IndexColorModel;
import java.awt.image.Raster;
import java.awt.image.WritableRaster;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;

public class PQCanvas extends Canvas {

	private static final long serialVersionUID = -5166271928949046848L;
	private PnnQuantizer pq = null;
	private Image image = null;

	public void set(final File file) {
		try {
			Image img = ImageIO.read(file);
			java.awt.MediaTracker tracker = new java.awt.MediaTracker(this);
			tracker.addImage(img, 0);
			try {
				tracker.waitForID(0);
			} catch (InterruptedException e) { }
			System.out.println("w = " + img.getWidth(this));
			System.out.println("h = " + img.getHeight(this));
			set(img);
		} catch (Exception ex) {
			ex.printStackTrace();
		}
	}
	
	private void addFileDrop() {
        new FileDrop(this, /*dragBorder,*/ new FileDrop.Listener() {
            public void filesDropped(java.io.File[] files) {
                try {
                    if (files.length > 0)
                        set(files[0]);
                } catch (Exception ex) {
                    java.util.logging.Logger.getLogger(PQCanvas.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
                    javax.swing.JOptionPane.showMessageDialog(PQCanvas.this, ex.getMessage(), "File not found", javax.swing.JOptionPane.ERROR_MESSAGE);
                }
            }   // end filesDropped
        }); // end FileDrop.Listener
    }
	
	private BufferedImage toIndexedBufferedImage(int[] pixels, IndexColorModel icm, int width, int height) {
	    //With this constructor we create an indexed bufferedimage with the same dimensiosn and with a default 256 color model
	    BufferedImage indexedImage= new BufferedImage(width, height,BufferedImage.TYPE_BYTE_INDEXED, icm);
	    WritableRaster raster = Raster.createWritableRaster(indexedImage.getSampleModel(), null);
		raster.setPixels(0, 0, width, height, pixels);
		indexedImage.setData(raster);
	    return indexedImage;
	}

	public void set(Image img) throws Exception {
		try {
			setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
			pq = new PnnLABQuantizer(img, this);
			int w = img.getWidth(this);
			int h = img.getHeight(this);
			int pixels[] = null;
			java.awt.image.PixelGrabber pg = new java.awt.image.PixelGrabber(img, 0, 0, w, h, pixels, 0, w);
			try {
				pg.grabPixels();
			} catch (InterruptedException e) { }
			if ((pg.getStatus() & java.awt.image.ImageObserver.ABORT) != 0)
				throw new IOException ("Image pixel grab aborted or errored");

			final int nMaxColors = 256;
			pixels = pq.convert(w, h, nMaxColors, true);
			if(pq.getColorModel() == null) {
				BufferedImage highColorImage= new BufferedImage(w, h, BufferedImage.TYPE_USHORT_565_RGB);
				short[] data = new short[pixels.length];
				for(int i=0; i<data.length; ++i)
					data[i] = (short) pixels[i];
				highColorImage.getRaster().setDataElements(0, 0, w, h, data);
			    this.image = highColorImage;
			}
			else
				this.image = toIndexedBufferedImage(pixels, pq.getColorModel(), w, h);
			//this.image = createImage(new MemoryImageSource(w, h, pixels, 0, w));
			this.getParent().setSize(w + 16, h + 38);
		} finally {
			setCursor(Cursor.getDefaultCursor());
		}
	}

	public void paint(Graphics graphics) {
		Graphics2D g2d = (Graphics2D) graphics.create();		
		if (image != null) {
			g2d.setRenderingHint(RenderingHints.KEY_INTERPOLATION, RenderingHints.VALUE_INTERPOLATION_BILINEAR); 
			graphics.drawImage(image, 0, 0, this);
		}
		else {
			graphics.setFont(new Font("Arial", Font.BOLD, 20));
			g2d.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_LCD_HRGB);
			graphics.drawString("Please drag an image file to here!", getWidth() / 5, getHeight() / 2);
		}
		g2d.dispose();
	} 

	public static void main(String [] args) throws java.io.IOException {
		PQCanvas canvas = new PQCanvas();
		java.awt.Frame frame = new java.awt.Frame("PnnQuant Test");
		frame.setSize(500, 500);
		canvas.addFileDrop();
		frame.add(canvas);
		Dimension dim = Toolkit.getDefaultToolkit().getScreenSize();
		frame.setLocation(dim.width/2-frame.getSize().width/2, dim.height/2-frame.getSize().height/2);
		frame.setVisible(true);
		frame.addWindowListener(new java.awt.event.WindowAdapter() {
			public void windowClosing(java.awt.event.WindowEvent e) { 
				System.out.println("Closing program");		
				System.exit(0); }
		} );

	}
}