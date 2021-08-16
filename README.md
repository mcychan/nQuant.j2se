# nQuant.j2se
<b>Fast pairwise nearest neighbor based algorithm with Java Swing</b>

Fast pairwise nearest neighbor based algorithm producing high quality 256 color 8 bit PNG images minimizing color loss for photo having red lips and supports 256 or less colors with transparency. 
A new multilevel thresholding algorithm derived from the well-known pairwise nearest neighbor ~PNN! method previously used in vector quantization.
The PNN has been considered a slow algorithm, as the time complexity of the original method was shown to be O(N3).<br/>
Faster approaches have been proposed by Kurita using a heap structure, and by Fra¨nti et al. using nearest neighbor pointers. The time complexity of the method is O(N log N), which is significantly better than that of optimal thresholding.

The Java version of the code, together with a test driver, using simple drag and drop method to convert the dropped image to 256 colors directly.
The main method is located at PQCanvas.java. It will open a JFrame dialog.

Original photo:<br><img src="https://i.stack.imgur.com/SE5x9.png" /><br>
Reduced to 256 colors by Pairwise Nearest Neighbor Quantization Algorithm:<br><img src="https://user-images.githubusercontent.com/26831069/129530423-0767e45a-001f-4024-aeed-49e897b5da62.png" />
