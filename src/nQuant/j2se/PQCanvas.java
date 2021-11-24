package nQuant.j2se;

import java.awt.BorderLayout;
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
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.DataFlavor;
import java.awt.event.ActionEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.net.URL;
import java.util.List;

import javax.imageio.ImageIO;
import javax.swing.AbstractAction;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.KeyStroke;
import javax.swing.SwingUtilities;
import javax.swing.SwingWorker;

public class PQCanvas extends Canvas {

	private static final long serialVersionUID = -5166271928949046848L;	
	private boolean hasAlpha;
	private BufferedImage image = null;
	private final TexturePaint tp;
	
	public PQCanvas() {
		tp = makeTexturePaint(16);
	}
	
	public void set(final BufferedImage img) throws Exception {
		if(img != null) {
			System.out.println("w = " + img.getWidth(this));
			System.out.println("h = " + img.getHeight(this));
			new PnnWorker(img).execute();
		}
		else
			javax.swing.JOptionPane.showMessageDialog(PQCanvas.this, "Cannot read image file", "Unknown image format", javax.swing.JOptionPane.ERROR_MESSAGE);
	}

	public void set(final File file) {
		try {
			set(ImageIO.read(file));			
		} catch (Exception ex) {
			ex.printStackTrace();
		}
	}
	
	public void set(final URL url) {
		try {
			set(ImageIO.read(url));
		} catch (Exception ex) {
			ex.printStackTrace();
		}
	}
	
	private void addClick() {
		addMouseListener(new MouseAdapter() {
			@Override
			public void mousePressed(MouseEvent e) {
				JFileChooser fileChooser = new JFileChooser();
				fileChooser.setCurrentDirectory(new File(System.getProperty("user.home")));
				if (fileChooser.showOpenDialog(PQCanvas.this) == JFileChooser.APPROVE_OPTION)
					set(fileChooser.getSelectedFile());
			}
		});
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

	private class PnnWorker extends SwingWorker<BufferedImage, String> { 
		private final BufferedImage img;

		PnnWorker(BufferedImage img) throws Exception {
			this.img = img;
		}

		@Override
		protected BufferedImage doInBackground() throws Exception {
			setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
			final PnnQuantizer pq = new PnnLABQuantizer(img, PQCanvas.this);
			BufferedImage highColorImage = pq.convert(256, true);
			hasAlpha = pq.hasAlpha();			
			return highColorImage;
			
			/* final Otsu otsu = new Otsu();
			BufferedImage bwImage = otsu.convertGrayScaleToBinary(img);
			hasAlpha = otsu.hasAlpha();			
			return bwImage; */
		}

		@Override
		protected void done() {
			try {
				image = get();								
			} catch (Exception e) {
				e.printStackTrace();
			} finally {	    	
				JFrame topFrame = (JFrame) SwingUtilities.getWindowAncestor(getParent());
				java.awt.Insets insets = topFrame.getInsets();
				topFrame.setSize(image.getWidth() + insets.right + insets.left, image.getHeight() + insets.top + insets.bottom);
				repaint();
				setCursor(Cursor.getDefaultCursor());
			}
		}
	};
	
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
		    java.awt.Insets insets = getParent().getInsets();
		    sz.setSize(sz.getWidth() - insets.right - insets.left, sz.getHeight() - insets.top - insets.bottom);
		    
		    if(hasAlpha) {       
		    	g2d.setPaint(tp);
		        g2d.fill(new Rectangle(sz));
		    }
		    		    
		    g2d.drawImage(image, null, 0, 0);
		}
		else {
			g2d.setFont(new Font("Arial", Font.BOLD, 20));
			g2d.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
			g2d.drawString("Please drag an image file to here!", getWidth() / 5, getHeight() / 2);
		}
		g2d.dispose();
	} 

	public static void main(String [] args) throws java.io.IOException {
		System.setProperty("java.net.useSystemProxies", "true");
		PQCanvas canvas = new PQCanvas();
		JFrame frame = new JFrame("PnnQuant Test");
		frame.setSize(500, 500);
		JPanel panel = new JPanel();
		panel.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW).put(KeyStroke.getKeyStroke("control V"), "paste");
        panel.getActionMap().put("paste", new AbstractAction() {
			private static final long serialVersionUID = 4914478601644487779L;

			public void actionPerformed(ActionEvent e) {
            	Clipboard clipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
                try {
                	DataFlavor[] keys = clipboard.getAvailableDataFlavors();
                	if(keys.length < 1)
                		return;
                	
                	Object data = clipboard.getData(keys[0]);
                    if (data instanceof String) {
                    	File file = new File((String) data);
                    	if(file.exists())
                    		canvas.set(file);
                    	else
                    		canvas.set(new URL((String) data));
                    }
                    else if(data instanceof BufferedImage)
                    	canvas.set((BufferedImage) data);
                    else if(data instanceof List)
                    	canvas.set(((List<File>) data).get(0));
				} catch (Exception ex) {
					java.util.logging.Logger.getLogger(PQCanvas.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
					javax.swing.JOptionPane.showMessageDialog(canvas, ex.getMessage(), "File not found", javax.swing.JOptionPane.ERROR_MESSAGE);
				}
            }
        });
		canvas.addClick();
		canvas.addFileDrop();
		panel.setLayout(new BorderLayout());
		panel.add(canvas, BorderLayout.CENTER);
		frame.setContentPane(panel);
		Dimension dim = Toolkit.getDefaultToolkit().getScreenSize();
		frame.setLocation(dim.width/2-frame.getSize().width/2, dim.height/2-frame.getSize().height/2);
		frame.setVisible(true);
		frame.addWindowListener(new java.awt.event.WindowAdapter() {
			public void windowClosing(java.awt.event.WindowEvent e) { 
				try {
					if(canvas.image != null) {
						String tempFilePath = System.getProperty("java.io.tmpdir") + "result.png";
						ImageIO.write(canvas.image, "png", new java.io.File(tempFilePath));
						System.out.println(tempFilePath);
					}
				} catch(Exception ex) {
					ex.printStackTrace();
				} finally {
					e.getWindow().dispose();
				}
			}
		});
	}
}