package nQuant.j2se;

import java.awt.Color;

public interface Ditherable {	
	public short nearestColorIndex(final Color[] palette, final Color c);
}
