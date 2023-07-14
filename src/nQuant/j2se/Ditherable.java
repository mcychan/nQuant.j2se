package nQuant.j2se;

import java.awt.Color;

public interface Ditherable {
	public int getColorIndex(final Color c);
	
	public short nearestColorIndex(final Color[] palette, final Color c, final int pos);
}
