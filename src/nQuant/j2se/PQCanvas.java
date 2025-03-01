package nQuant.j2se;

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
import java.awt.datatransfer.Transferable;
import java.awt.event.ActionEvent;
import java.awt.event.InputEvent;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.awt.geom.AffineTransform;
import java.awt.geom.Rectangle2D;
import java.awt.image.AffineTransformOp;
import java.awt.image.BufferedImage;
import java.awt.image.IndexColorModel;
import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.net.URL;
import java.security.SecureRandom;
import java.security.cert.CertificateException;
import java.security.cert.X509Certificate;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import javax.imageio.ImageIO;
import javax.net.ssl.HostnameVerifier;
import javax.net.ssl.HttpsURLConnection;
import javax.net.ssl.SSLContext;
import javax.net.ssl.SSLSession;
import javax.net.ssl.X509TrustManager;
import javax.swing.AbstractAction;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JViewport;
import javax.swing.KeyStroke;
import javax.swing.Scrollable;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import javax.swing.SwingWorker;
import javax.swing.TransferHandler;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import nQuant.ga.APNsgaIII;
import nQuant.ga.NsgaIII;

public class PQCanvas extends JPanel implements Scrollable, MouseWheelListener {

	private static final long serialVersionUID = -5166271928949046848L;	
	private volatile boolean hasAlpha;
	private int maxUnitIncrement = 1;
	private double zoom = 1.0;
	private BufferedImage image = null;
	private List<BufferedImage> images = null;
	private final TexturePaint tp;
	
	public PQCanvas() {
		tp = makeTexturePaint(16);
		trustEveryone();
	}
	
	public void set(final List<BufferedImage> imgs) throws Exception {
		if(!imgs.isEmpty())	
			new PnnWorker(imgs).execute();
		else
			javax.swing.JOptionPane.showMessageDialog(PQCanvas.this, "Cannot read image", "Unknown image format", javax.swing.JOptionPane.ERROR_MESSAGE);
	}

	public void set(final File[] files) {
		try {
			List<BufferedImage> images = new ArrayList<>();
			for(File file : files) {
				BufferedImage img = ImageIO.read(file);
				if(img == null)
					throw new UnsupportedOperationException(file.getAbsolutePath());
				images.add(img);
			}
			set(images);
		} catch (Exception ex) {
			ex.printStackTrace();
		}
	}
	
	public void set(final URL url) {
		try {
			BufferedImage img = ImageIO.read(url);
			if(img == null)
				throw new UnsupportedOperationException(url.getPath());
			set(Collections.singletonList(img));
		} catch (Exception ex) {
			ex.printStackTrace();
		}
	}	
	
	public double getZoom() {
		return zoom;
	}

	public void setZoom(double zoom) {
		if(zoom > .05 && zoom < 100)
			this.zoom = zoom;
	}
	
	public Dimension getPreferredSize() {
		if(image != null)
			return new Dimension((int) Math.ceil(image.getWidth() * getZoom()), (int) Math.ceil(image.getHeight() * getZoom()));
		return getParent().getSize();
	}

	public Dimension getPreferredScrollableViewportSize() {
		return getPreferredSize();
	}

	public int getScrollableUnitIncrement(Rectangle visibleRect, int orientation, int direction) {
		//Get the current position.
		int currentPosition = (orientation == SwingConstants.HORIZONTAL) ? visibleRect.x : visibleRect.y;

		//Return the number of pixels between currentPosition
		//and the nearest tick mark in the indicated direction.
		if (direction < 0) {
			int newPosition = currentPosition - (currentPosition / maxUnitIncrement) * maxUnitIncrement;
			return (newPosition == 0) ? maxUnitIncrement : newPosition;
		}
		
		return ((currentPosition / maxUnitIncrement) + 1) * maxUnitIncrement - currentPosition;
	}

	public int getScrollableBlockIncrement(Rectangle visibleRect, int orientation, int direction) {
		if (orientation == SwingConstants.HORIZONTAL)
			return visibleRect.width - maxUnitIncrement;
		return visibleRect.height - maxUnitIncrement;
	}

	public boolean getScrollableTracksViewportWidth() {
		return false;
	}

	public boolean getScrollableTracksViewportHeight() {
		return false;
	}

	public void setMaxUnitIncrement(int pixels) {
		maxUnitIncrement = pixels;
	}

	private void addClick() {
		addMouseListener(new MouseAdapter() {
			@Override
			public void mouseReleased(MouseEvent e) {
				JFileChooser fileChooser = new JFileChooser();
				fileChooser.setCurrentDirectory(new File(System.getProperty("user.home")));
				if (fileChooser.showOpenDialog(PQCanvas.this) == JFileChooser.APPROVE_OPTION)
					set(fileChooser.getSelectedFiles());
			}
		});
	}

	@Override
	public TransferHandler getTransferHandler() {
		return new TransferHandler() {
			@Override
			public boolean canImport(JComponent component, DataFlavor[] flavors) {
				return true;
			}
			
			@Override
			public boolean importData(JComponent component, Transferable transferable) {
				try {
					for (DataFlavor flavor : transferable.getTransferDataFlavors()) {						
						if (DataFlavor.imageFlavor.equals(flavor)) {
							BufferedImage img = (BufferedImage) transferable.getTransferData(DataFlavor.imageFlavor);
							set(Collections.singletonList(img));
							return true;
						}
						List<File> files = (List<File>) transferable.getTransferData(DataFlavor.javaFileListFlavor);
						if (files != null) {
							if (files.size() > 0) {
								set(files.toArray(new File[0]));
								return true;
							}
						}
						if (DataFlavor.stringFlavor.equals(flavor)) {
							URL url = new URL((String) transferable.getTransferData(DataFlavor.stringFlavor));
							set(url);
							return true;
						}
					}
				} catch (Exception ex) {
					java.util.logging.Logger.getLogger(PQCanvas.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
					javax.swing.JOptionPane.showMessageDialog(PQCanvas.this, ex.getMessage(), "File not found", javax.swing.JOptionPane.ERROR_MESSAGE);
				}
				return false;
			}
			
			@Override
			public int getSourceActions(JComponent component) {
				return COPY;
			}
		};
	}

	private class PnnWorker extends SwingWorker<BufferedImage, String> { 
		private final List<BufferedImage> imgs;
		private long startTime;

		PnnWorker(List<BufferedImage> imgs) throws Exception {
			this.imgs = imgs;
		}

		@Override
		protected BufferedImage doInBackground() throws Exception {
			setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
			startTime = System.currentTimeMillis();
			final PnnQuantizer pq = new PnnLABQuantizer(imgs.get(0), PQCanvas.this);
			BufferedImage highColorImage = pq.convert(256, true);
			hasAlpha = pq.hasAlpha();
			return highColorImage;

			/* NsgaIII<PnnLABGAQuantizer> alg = new APNsgaIII<>(new PnnLABGAQuantizer(new PnnLABQuantizer(imgs.get(0), PQCanvas.this), imgs, 256), 2, 2, 80, 3);
			alg.run(9999, -Math.ulp(1.0));
			try (PnnLABGAQuantizer pGAq = alg.getResult()) {
				System.out.println("\n" + pGAq.getResult());
				images = pGAq.convert(true);
				hasAlpha = pGAq.hasAlpha();
				return images.get(0);
			} */

			/* final Otsu otsu = new Otsu();
			BufferedImage bwImage = otsu.convertGrayScaleToBinary(imgs.get(0));
			hasAlpha = otsu.hasAlpha();
			return bwImage; */
		}

		@Override
		protected void done() {
			try {
				image = get();
				System.out.println("\n" + (System.currentTimeMillis() - startTime) / 1000.0 + " seconds");
			} catch (Exception e) {
				e.printStackTrace();
			} finally {
				final JFrame topFrame = (JFrame) SwingUtilities.getWindowAncestor(getParent());
				topFrame.setVisible(false);
				java.awt.Insets insets = topFrame.getInsets();
				final Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
				if(image.getWidth() > screenSize.getWidth() || image.getHeight() > screenSize.getHeight()) {					
					topFrame.setLocation(-insets.left, 0);
					topFrame.setSize((int) screenSize.getWidth() + insets.right + insets.left, (int) screenSize.getHeight() - insets.top - insets.bottom);
					setZoom(Math.min(screenSize.getWidth() / image.getWidth(), screenSize.getHeight() / image.getHeight()));
				}
				else {
					topFrame.setSize((int) (image.getWidth() + insets.right + 1.5 * insets.left), (int) (image.getHeight() + insets.top + 1.5 * insets.bottom));
					setZoom(Math.min(topFrame.getWidth() / image.getWidth(), topFrame.getHeight() / image.getHeight()));
				}

				final JScrollPane scrollPane = (JScrollPane) SwingUtilities.getAncestorOfClass(JScrollPane.class, getParent());
				scrollPane.getHorizontalScrollBar().setValue(0);
				scrollPane.getVerticalScrollBar().setValue(0);
				repaint();
				topFrame.setVisible(true);
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

	@Override
	public void paintComponent(Graphics g) {
		super.paintComponent(g);
		Graphics2D g2d = (Graphics2D) g.create();
		g2d.setRenderingHint(RenderingHints.KEY_INTERPOLATION, RenderingHints.VALUE_INTERPOLATION_BICUBIC);		

		if (image != null) {
			g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_OFF);		    
			g2d.setRenderingHint(RenderingHints.KEY_ALPHA_INTERPOLATION, RenderingHints.VALUE_ALPHA_INTERPOLATION_QUALITY);
			g2d.setRenderingHint(RenderingHints.KEY_COLOR_RENDERING, RenderingHints.VALUE_COLOR_RENDER_DEFAULT);
			g2d.setRenderingHint(RenderingHints.KEY_DITHERING, RenderingHints.VALUE_DITHER_DISABLE);
			g2d.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_SPEED);
			
			if(hasAlpha) {
				g2d.setPaint(tp);
				g2d.fill(new Rectangle(getSize()));
			}

		    AffineTransform scaleInstance = AffineTransform.getScaleInstance(zoom, zoom);
		    AffineTransformOp scaleOp = new AffineTransformOp(scaleInstance, AffineTransformOp.TYPE_BILINEAR);		    

			g2d.drawImage(image, scaleOp, 0, 0);
		}
		else {
			g2d.setFont(new Font("Arial", Font.BOLD, 20));
			g2d.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
			g2d.clearRect(0, 0, getParent().getWidth(), getParent().getHeight());
			g2d.drawString("Please drag an image file to here!", getWidth() / 5, getHeight() / 2);
		}
		g2d.dispose();
	}
	
	@Override
	public void mouseWheelMoved(MouseWheelEvent e) {
		JScrollPane scrollPane = (JScrollPane) e.getSource();
		if ((e.isControlDown() || !scrollPane.getVerticalScrollBar().isVisible()) && image != null) {
			double oldZoom = getZoom();
			double amount = Math.pow(1.1, e.getScrollAmount());
			if (e.getWheelRotation() > 0) {
				//zoom out (amount)
				setZoom(oldZoom / amount);
			} else {
				//zoom in (amount)
				setZoom(oldZoom * amount);
			}
			
			if(Math.abs(oldZoom - getZoom()) > 0.05) {
				scrollPane.getViewport().setSize((int)(image.getWidth() * getZoom()), (int)(image.getHeight() * getZoom()));
				final JFrame topFrame = (JFrame) SwingUtilities.getWindowAncestor(scrollPane);
				java.awt.Insets insets = topFrame.getInsets();
				final Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
				final Dimension scrollSize = scrollPane.getViewport().getSize();
				if(scrollSize.getWidth() > screenSize.getWidth() || scrollSize.getHeight() > screenSize.getHeight()) {
					topFrame.setSize((int) screenSize.getWidth() + insets.right + insets.left, (int) screenSize.getHeight() - insets.top - insets.bottom);
					topFrame.setLocation(-insets.left, 0);
				}
				else
					topFrame.setSize((int) (scrollSize.getWidth() + insets.right + 1.5 * insets.left), (int) (scrollSize.getHeight() + insets.top + 1.5 * insets.bottom));

				repaint();
			}
		}
		else {
			// if ctrl isn't down then propagate event to parent
			scrollPane.getParent().dispatchEvent(e);
		}
	}

	private static void trustEveryone() { 
		try { 
			HttpsURLConnection.setDefaultHostnameVerifier(new HostnameVerifier() {
				public boolean verify(String hostname, SSLSession session) {
					return true; 
				}}); 
			SSLContext context = SSLContext.getInstance("TLS"); 
			context.init(null, new X509TrustManager[]{new X509TrustManager() {
				public void checkClientTrusted(X509Certificate[] chain, String authType) throws CertificateException {
					if(chain.length < 1)
						throw new CertificateException();
				}
				public void checkServerTrusted(X509Certificate[] chain, String authType) throws CertificateException {
					if(chain.length < 1)
						throw new CertificateException();
				}
				public X509Certificate[] getAcceptedIssuers() {
					return new X509Certificate[0];
				}}}, new SecureRandom()); 
			HttpsURLConnection.setDefaultSSLSocketFactory(context.getSocketFactory());
		} catch (Exception e) { // should never happen 
			e.printStackTrace(); 
		}
	}

	public static void main(String [] args) throws java.io.IOException {
		System.setProperty("java.net.useSystemProxies", "true");
		PQCanvas canvas = new PQCanvas();
		JFrame frame = new JFrame("PnnQuant Playground");
		frame.setPreferredSize(new Dimension(500, 500));
		JScrollPane scrollPane = new JScrollPane(canvas);
		scrollPane.getViewport().addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e){
				JViewport viewport = (JViewport) e.getSource();
				viewport.revalidate();
			}
		});
		scrollPane.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW).put(KeyStroke.getKeyStroke(KeyEvent.VK_V, InputEvent.CTRL_DOWN_MASK), "paste");
		scrollPane.getActionMap().put("paste", new AbstractAction() {
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
							canvas.set(new File[] {file});
						else
							canvas.set(new URL((String) data));
					}
					else if(data instanceof BufferedImage)
						canvas.set(Collections.singletonList((BufferedImage) data));
					else if(data instanceof List)
						canvas.set(((List<File>) data).toArray(new File[0]));
				} catch (Exception ex) {
					java.util.logging.Logger.getLogger(PQCanvas.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
					javax.swing.JOptionPane.showMessageDialog(canvas, ex.getMessage(), "File not found", javax.swing.JOptionPane.ERROR_MESSAGE);
				}
			}
		});
		canvas.addClick();
		canvas.setTransferHandler(canvas.getTransferHandler());
		
		scrollPane.addMouseWheelListener(canvas);
		scrollPane.getViewport().setScrollMode(JViewport.SIMPLE_SCROLL_MODE);
		frame.pack();
		frame.setContentPane(scrollPane);
		frame.setLocationRelativeTo(null);
		frame.setVisible(true);
		frame.addWindowListener(new java.awt.event.WindowAdapter() {
			public void windowClosing(java.awt.event.WindowEvent e) { 
				try {
					if(canvas.image != null) {
						final String fileType = canvas.images != null || canvas.image.getColorModel() instanceof IndexColorModel ? "gif" : "png";
						String tempFilePath = System.getProperty("java.io.tmpdir") + "result." + fileType;
						if(canvas.images != null) {
							try (GifWriter writer = new GifWriter(tempFilePath, canvas.hasAlpha, 850, true)) {
								canvas.images.stream().forEach(img -> {
									try {
										writer.writeToSequence(img);
									} catch (IOException e1) {
										throw new UncheckedIOException(e1);
									}
								});
							}
						}
						else
							ImageIO.write(canvas.image, fileType, new java.io.File(tempFilePath));
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