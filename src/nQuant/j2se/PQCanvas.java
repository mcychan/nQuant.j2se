package nQuant.j2se;

import java.awt.Canvas;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.RenderingHints;
import java.awt.TexturePaint;
import java.awt.Toolkit;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.awt.image.ColorModel;
import java.awt.image.DataBufferShort;
import java.awt.image.DirectColorModel;
import java.awt.image.IndexColorModel;
import java.awt.image.Raster;
import java.awt.image.RenderedImage;
import java.awt.image.WritableRaster;
import java.io.File;

import javax.imageio.ImageIO;
import javax.swing.SwingWorker;

public class PQCanvas extends Canvas {

	private static final long serialVersionUID = -5166271928949046848L;	
	private boolean hasAlpha;
	private BufferedImage image = null;
	private final TexturePaint tp;
	
	public PQCanvas() {
		tp = makeTexturePaint(16);
	}

	public void set(final File file) {
		try {
			BufferedImage img = ImageIO.read(file);
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
		WritableRaster raster = Raster.createWritableRaster(icm.createCompatibleSampleModel(width, height), new DataBufferShort(qPixels, qPixels.length), null);
		return new BufferedImage(icm, raster, icm.isAlphaPremultiplied(), null);
	}

	private class PnnWorker extends SwingWorker<BufferedImage, String> { 
		private final PnnQuantizer pq;

		PnnWorker(PnnQuantizer pq) {
			this.pq = pq;
		}

		@Override
		protected BufferedImage doInBackground() throws Exception {
			setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
			int w = pq.getWidth();
			int h = pq.getHeight();

			short[] qPixels = pq.convert(256, true);
			if(pq.getColorModel() instanceof IndexColorModel)
				return toIndexedBufferedImage(qPixels, (IndexColorModel) pq.getColorModel(), w, h);

			BufferedImage highColorImage = null;
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
				hasAlpha = pq.hasAlpha();
				String tempFilePath = System.getProperty("java.io.tmpdir") + "result.png";
				ImageIO.write((RenderedImage) image, "png", new java.io.File(tempFilePath));
				System.out.println(tempFilePath);
			} catch (Exception e) {
				e.printStackTrace();
			} finally {	    	
				getParent().setSize(pq.getWidth() + 16, pq.getHeight() + 38);
				repaint();
				setCursor(Cursor.getDefaultCursor());
			}
		}
	};

	public void set(BufferedImage img) throws Exception {
		PnnQuantizer pq = new PnnLABQuantizer(img, this);
		new PnnWorker(pq).execute();		
	}
	
	private TexturePaint makeTexturePaint(int size) {
		BufferedImage img = new BufferedImage(size, size, BufferedImage.TYPE_INT_ARGB);
		java.awt.Color[] colors = new java.awt.Color[] {new java.awt.Color(255, 255, 255, 127), new java.awt.Color(192, 192, 192, 127)};
        
		int s2 = size / 2;
		Graphics2D g2d = img.createGraphics();	        
	    g2d.setColor(colors[0]);
	    g2d.fillRect(0, 0, s2, s2);
	    g2d.setColor(colors[1]);
	    g2d.fillRect(0, s2, s2, s2);
	    g2d.fillRect(s2, 0, s2, s2);
	    g2d.setColor(colors[0]);
	    g2d.fillRect(s2, s2, s2, s2);
	    g2d.dispose();
	    Rectangle2D bounds = new Rectangle2D.Float(0, 0, size, size);
	    return new TexturePaint(img, bounds);
	}

	public void paint(Graphics graphics) {
		Graphics2D g2d = (Graphics2D) graphics;
		g2d.setRenderingHint(RenderingHints.KEY_INTERPOLATION, RenderingHints.VALUE_INTERPOLATION_BICUBIC);
		if (image != null) {			
			g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_OFF);		    
		    g2d.setRenderingHint(RenderingHints.KEY_ALPHA_INTERPOLATION, RenderingHints.VALUE_ALPHA_INTERPOLATION_QUALITY);
		    g2d.setRenderingHint(RenderingHints.KEY_COLOR_RENDERING, RenderingHints.VALUE_COLOR_RENDER_DEFAULT);
		    g2d.setRenderingHint(RenderingHints.KEY_DITHERING, RenderingHints.VALUE_DITHER_DISABLE);
		    g2d.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_SPEED);
		    
		    Dimension sz = getParent().getSize();
		    sz.setSize(sz.getWidth() - 16, sz.getHeight() - 38);
		    
		    if(hasAlpha) {       
		    	g2d.setPaint(tp);
		        g2d.fill(new Rectangle(sz));
		    }
		    		    
		    g2d.drawImage(image, 0, 0, (int) sz.getWidth(), (int) sz.getHeight(), this);
		}
		else {
			g2d.setFont(new Font("Arial", Font.BOLD, 20));
			g2d.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
			g2d.drawString("Please drag an image file to here!", getWidth() / 5, getHeight() / 2);
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
		});
	}
}