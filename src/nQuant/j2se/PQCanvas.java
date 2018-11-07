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
import java.awt.image.ColorModel;
import java.awt.image.DataBufferByte;
import java.awt.image.DirectColorModel;
import java.awt.image.IndexColorModel;
import java.awt.image.Raster;
import java.awt.image.WritableRaster;
import java.io.File;

import javax.imageio.ImageIO;
import javax.swing.SwingWorker;

public class PQCanvas extends Canvas {

	private static final long serialVersionUID = -5166271928949046848L;
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

	private BufferedImage toIndexedBufferedImage(short[] qPixels, IndexColorModel icm, int width, int height) {
		//With this constructor we create an indexed bufferedimage with the same dimensiosn and with a default 256 color model
		BufferedImage indexedImage= new BufferedImage(width, height,BufferedImage.TYPE_BYTE_INDEXED, icm);
		byte[] data = new byte[qPixels.length];
		for(int i=0; i<data.length; ++i)
			data[i] = (byte) qPixels[i];
		WritableRaster raster = Raster.createWritableRaster(indexedImage.getSampleModel(), new DataBufferByte(data, data.length), null);
		indexedImage.setData(raster);
		return indexedImage;
	}

	private class PnnWorker extends SwingWorker<Image, String> { 
		private final PnnQuantizer pq;

		PnnWorker(PnnQuantizer pq) {
			this.pq = pq;
		}

		@Override
		protected Image doInBackground() throws Exception {
			setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
			int w = pq.getWidth();
			int h = pq.getHeight();

			short[] qPixels = pq.convert(w, h, 256, true);
			if(pq.getColorModel() instanceof IndexColorModel)
				return toIndexedBufferedImage(qPixels, (IndexColorModel) pq.getColorModel(), w, h);

			BufferedImage highColorImage= null;
			if(pq.getColorModel() instanceof DirectColorModel) {
				ColorModel cmSw = pq.getColorModel();
				WritableRaster wr = cmSw.createCompatibleWritableRaster(w, h);
				highColorImage = new BufferedImage(cmSw, wr, cmSw.isAlphaPremultiplied(), null);
			}
			else
				highColorImage = new BufferedImage(w, h, BufferedImage.TYPE_USHORT_565_RGB);

			highColorImage.getRaster().setDataElements(0, 0, w, h, qPixels);
			return highColorImage;
		}

		@Override
		protected void done() {
			try {
				image = get();
			} catch (Exception e) {
				e.printStackTrace();
			} finally {	    	
				getParent().setSize(pq.getWidth() + 16, pq.getHeight() + 38);
				repaint();
				setCursor(Cursor.getDefaultCursor());
			}
		}
	};

	public void set(Image img) throws Exception {
		PnnQuantizer pq = new PnnLABQuantizer(img, this);
		new PnnWorker(pq).execute();		
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