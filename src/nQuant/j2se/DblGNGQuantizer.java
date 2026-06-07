package nQuant.j2se;

/* Distributed Batch Learning Growing Neural Gas algorithm with CIELAB color space
Copyright (c) 2026 Miller Cy Chan
 * Siow, C. Z., Saputra, A. A., Obo, T., & Kubota, N. (2024).
 * Distributed batch learning of growing neural gas for quick and efficient clustering.
 * Mathematics, 12(12), 1909. */

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.awt.image.ColorModel;
import java.awt.image.ImageObserver;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.stream.Collectors;

import nQuant.j2se.CIELABConvertor.Lab;
import nQuant.j2se.CIELABConvertor.MutableDouble;

import static nQuant.j2se.BitmapUtilities.BYTE_MAX;

public class DblGNGQuantizer {
	private byte alphaThreshold = 0xF;
	private int maxNodes = 256; // max nodes per model
	private int epochs = 20, maxAge = 10; // Maximum age in terms of epochs/batches, not pixels

	private final List<GNGNode> nodes = new ArrayList<>();
	private static Random random;

	protected boolean enforcedDither = true, hasSemiTransparency = false;
	protected int m_transparentPixelIndex = -1;
	private static final double TRANS_RATE = 1 - (512 + 101) / 768.0;
	protected final int width, height;
	protected int[] pixels;
	protected Color m_transparentColor = new Color(BYTE_MAX, BYTE_MAX, BYTE_MAX, 0);
	protected double PR = 0.299, PG = 0.587, PB = 0.114, PA = .3333;
	protected static final float[][] coeffs = new float[][] {
		{0.299f, 0.587f, 0.114f},
		{-0.14713f, -0.28886f, 0.436f},
		{0.615f, -0.51499f, -0.10001f}
	};
	protected float[] saliencies;

	protected Color[] m_palette;
	private ColorModel m_colorModel;
	protected Map<Integer, Lab> pixelMap = new HashMap<>();
	protected Map<Integer, int[]> closestMap = new HashMap<>();
	protected Map<Integer, Short> nearestMap = new HashMap<>();

	private DblGNGQuantizer(BufferedImage im, int w, int h) {
		width = w;
		height = h;
		setPixels(im);
	}

	public DblGNGQuantizer(BufferedImage im, ImageObserver obs) {
		this(im, im.getWidth(obs), im.getHeight(obs));
	}

	private void setPixels(BufferedImage im) {
		pixels = im.getRGB(0, 0, width, height, null, 0, width);
	}
	
	protected int getColorIndex(final Color c)
	{
		if(hasSemiTransparency)
			return (c.getAlpha() & 0xF0) << 8 | (c.getRed() & 0xF0) << 4 | (c.getGreen() & 0xF0) | (c.getBlue() >> 4);
		if(m_transparentPixelIndex > -1)
			return (c.getAlpha() & 0x80) << 8 | (c.getRed() & 0xF8) << 7 | (c.getGreen() & 0xF8) << 2 | (c.getBlue() >> 3);
		return (c.getRed() & 0xF8) << 8 | (c.getGreen() & 0xFC) << 3 | (c.getBlue() >> 3);
	}

	Lab getLab(final int pixel)
	{
		Lab lab1 = pixelMap.get(pixel);
		if (lab1 == null) {
			lab1 = CIELABConvertor.RGB2LAB(pixel, hasAlpha());
			pixelMap.put(pixel, lab1);
		}
		return lab1;
	}

	protected final void setColorModel(final Color[] palette)
	{
		m_palette = palette;
		m_colorModel = BitmapUtilities.createColorModel(palette, m_transparentPixelIndex, hasSemiTransparency, hasSemiTransparency);
	}

	static class GNGNode {
		final double[] weight;	// LAB(alpha) vector
		double error = 0.0;
		private final Map<GNGNode, Integer> neighbors = new HashMap<>();

		public GNGNode(double[] weight) {
			this.weight = weight;		
		}
		
		public void addNeighbour(GNGNode nextNode) {
			neighbors.put(nextNode, 0);
		}
		
		public GNGNode findNeighborByMaxError() {
			return neighbors.keySet().stream().max(Comparator.comparingDouble(n -> n.error))
			.orElse(null);
		}
		
		public void incrementAge() {
			neighbors.replaceAll((neighbor, age) -> age + 1);
		}
		
		public boolean noNeighbor() {
			return neighbors.isEmpty();
		}
		
		public void removeNeighbour(final GNGNode nextNode) {
			neighbors.keySet().remove(nextNode);
		}
		
		public void removeNeighbourByAge(final int maxAge) {
			neighbors.values().removeIf(age -> age > maxAge);
		}

		public double distance(double[] input) {
			double d = 0;
			for (int i = 0; i < 3; ++i) {
				double diff = weight[i] - input[i];
				d += diff * diff;
			}

			if (weight.length > 3) {
				d += TRANS_RATE * BitmapUtilities.sqr(weight[3] - input[3]);
			}
			return d;
		}
	}

	private int calculateStartingPoints(double learningRate) {
		// k controls the steepness of the curve. 10.0 is optimal for 32 maxNodes.
		final double K = hasAlpha() || maxNodes > 32 ? 6.5 : 10.0; 

		// Calculate the exponential relationship
		double continuousPoints = maxNodes * Math.exp(-K * learningRate);

		int noOfStartingPoints = (int) Math.round(continuousPoints);

		// Guardrails: 
		// Absolute floor must be 2 (GNG requires a minimum of 1 edge / 2 nodes to start).
		// Absolute ceiling should be capped at 25% of maxNodes so the network has room to "Grow".
		// 4D TRANSPARENCY GUARDRAILS
	    // - Floor (8): Ensures the network can simultaneously track solid background colors, 
	    //   foreground assets, and alpha gradient transitions from epoch zero.
		//   leaving 75% of the node budget open for organic GNG growth.
		int minFloor = hasAlpha() ? 8 : 2;
		int maxCeiling = Math.max(2, maxNodes / 4); 

		return Math.max(minFloor, Math.min(maxCeiling, noOfStartingPoints));
	}

	public void initializeDistributedNode(List<GNGNode> samples, int noOfStartingPoints) {
		if (samples == null || samples.isEmpty()) {
			throw new IllegalArgumentException("Sample list cannot be empty.");
		}

		// Safety check: Cannot initialize more starting points than available pixels
		int actualStartingPoints = Math.min(noOfStartingPoints, samples.size());
		// GNG requires at least 2 nodes to form initial structural edges
		if (actualStartingPoints < 2) {
			actualStartingPoints = 2;
		}

		maxAge = Math.max(maxAge, actualStartingPoints * (epochs / 10));

		random = new Random(samples.size());
		nodes.clear();

		// 1. Collect random unique initial points from sample pool
		List<GNGNode> initialNodes = new ArrayList<>();
		Set<Integer> chosenIndices = new HashSet<>();

		while (initialNodes.size() < actualStartingPoints) {
			int randomIndex = random.nextInt(samples.size());
			// Ensure we don't pick the exact same pixel index twice
			if (chosenIndices.add(randomIndex)) {
				initialNodes.add(samples.get(randomIndex));
			}
		}

		// 2. Link initial nodes linearly to establish the starting graph topology
		// GNG relies on neighbor relationships to start learning and prevent pruning
		for (int i = 0; i < initialNodes.size(); ++i) {
			GNGNode currentNode = initialNodes.get(i);

			// Link to the next node in the list (wrap around at the end)
			GNGNode nextNode = initialNodes.get((i + 1) % initialNodes.size());

			currentNode.addNeighbour(nextNode);
			nextNode.addNeighbour(currentNode);
		}

		// 3. Batch transfer into the active thread-safe node list
		nodes.addAll(initialNodes);
	}

	private void insertNewNodeWeighted(Map<GNGNode, List<GNGNode>> assignments) {
		GNGNode q = nodes.stream()
			.max(Comparator.comparingDouble(n -> {
				List<GNGNode> samples = assignments.get(n);
				if (samples == null || samples.isEmpty())
					return 0.0;
	
				// Optimize for average error per pixel rather than total cluster error
				return n.error / Math.log1p(samples.size()); 
			}))
			.orElse(null);

		if (q == null || q.noNeighbor())
			return;

		GNGNode f = q.findNeighborByMaxError();

		if (f == null)
			return;

		double[] newWeight = new double[q.weight.length];
		for (int i = 0; i < q.weight.length; ++i) {
			newWeight[i] = (q.weight[i] + f.weight[i]) / 2.0;
		}

		GNGNode r = new GNGNode(newWeight);

		q.removeNeighbour(f);
		f.removeNeighbour(q);

		q.addNeighbour(r);
		f.addNeighbour(r);
		r.addNeighbour(q);
		r.addNeighbour(f);

		nodes.add(r);

		q.error *= 0.5;
		f.error *= 0.5;
		r.error = q.error;
	}

	/**
	 * Adapts node weights based on a progressive learning schedule.
	 * * @param assignments The map containing samples grouped under their winning GNGNode.
	 * @param baseLearningRate The dynamic learning rate calculated for the epoch loop.
	 * @param progress The current epoch progress expressed as a fraction from 0.0 to 1.0.
	 */
	private void updateNodeWeightsAdaptive(
		Map<GNGNode, List<GNGNode>> assignments, 
		double baseLearningRate, 
		double progress
	) {
		// Process each node cluster in parallel for high performance
		assignments.entrySet().parallelStream().forEach(entry -> {
			GNGNode node = entry.getKey();
			List<GNGNode> cluster = entry.getValue();
			
			// Safety check: skip if no samples were mapped to this node during the batch
			if (cluster == null || cluster.isEmpty())
				return;

			// 1. Calculate the arithmetic mean (centroid) of the assigned color samples
			double[] mean = new double[node.weight.length];
			for (GNGNode sample : cluster) {
				for (int i = 0; i < mean.length; ++i) {
					mean[i] += sample.weight[i];
				}
			}
			for (int i = 0; i < mean.length; ++i) {
				mean[i] /= cluster.size();
			}

			// 2. Apply the Adaptive Schedule
			if (progress < 0.4) {
				// WARM-UP EXPLORATION: Decay the learning rate slightly over time 
				// so nodes glide smoothly into position, preventing catastrophic hue shifts.
				double currentLR = baseLearningRate * (1.0 - progress); 
				
				for (int i = 0; i < node.weight.length; ++i) {
					node.weight[i] += currentLR * (mean[i] - node.weight[i]);
				}
			} else {
				// FINE-TUNING EXPLOITATION: Linearly scale up the snapping factor from 0.0 to 1.0.
				// This pulls the nodes precisely onto the true cluster centers, eliminating color bleeding.
				double snapFactor = (progress - 0.4) / 0.6; 
				
				for (int i = 0; i < node.weight.length; ++i) {
					// Linear interpolation (LERP) between current weight and true mean
					node.weight[i] = node.weight[i] + snapFactor * (mean[i] - node.weight[i]);
				}
			}
		});
	}
	
	/**
	 * Manages edge ages, prunes dead connections, and handles deterministic node insertion.
	 * * @param assignments The map containing samples grouped under their winning GNGNode.
	 * @param maxNodes The target ceiling for your color palette size (e.g., 32).
	 * @param remainingEpochs How many epochs are left in the growth phase (used to scale insertion).
	 */
	private void manageGraphTopology(
		Map<GNGNode, List<GNGNode>> assignments, 
		int remainingEpochs
	) {
		// Hyperparameter configuration
		final int MAX_AGE = 20; // Maximum allowed age for an edge before it is pruned

		// 1. AGE INCREMENTATION
		// Increment the age of all topological connections currently in the graph
		for (GNGNode node : nodes) {
			node.incrementAge();
		}

		// 2. RECONNECT ACTIVE EDGES
		// For every node that claimed samples this epoch, make sure its connection 
		// to the runner-up node is refreshed back to age 0.
		assignments.forEach((firstWinner, samplesInCluster) -> {
			if (samplesInCluster == null || samplesInCluster.isEmpty())
				return;

			// Find the absolute closest sample in this bucket to act as our topological anchor
			double[] anchorSample = samplesInCluster.get(0).weight; 

			// Find the second closest node in the entire network to this anchor sample
			GNGNode secondWinner = null;
			double minDistance2 = Double.MAX_VALUE;

			for (GNGNode potentialSecond : nodes) {
				if (potentialSecond == firstWinner) continue; // Cannot pair with itself
				
				double dist = potentialSecond.distance(anchorSample);
				if (dist < minDistance2) {
					minDistance2 = dist;
					secondWinner = potentialSecond;
				}
			}

			// Refresh or create the link between the two closest matching colors
			if (secondWinner != null) {
				firstWinner.addNeighbour(secondWinner);
				secondWinner.addNeighbour(firstWinner);
			}
		});

		// 3. PRUNE DEAD CONNECTIONS
		// Wipe out edges that haven't been validated recently
		for (GNGNode node : nodes) {
			node.removeNeighbourByAge(MAX_AGE);
		}

		// 4. CRITICAL RESCUE CHECK
		// Only remove a node if it has lost all neighbors AND failed to capture 
		// a single sample point. This preserves fragile, low-frequency minor color accents.
		nodes.removeIf(node -> {
			if (node.noNeighbor()) {
				List<GNGNode> assigned = assignments.get(node);
				return assigned == null || assigned.isEmpty();
			}
			return false;
		});

		// 5. DETERMINISTIC NODE GROWTH
		// Dynamically calculate how many nodes to insert to steadily reach maxNodes 
		// by the end of Phase 2, avoiding a sudden rush of insertions at the end.
		int missingNodes = maxNodes - nodes.size();
		if (missingNodes > 0 && remainingEpochs > 0) {
			// Linearly distribute the required insertions over the remaining growth epochs
			int targetInsertions = (int) Math.ceil((double) missingNodes / remainingEpochs);
			
			for (int i = 0; i < targetInsertions; ++i) {
				if (nodes.size() < maxNodes) {
					insertNewNodeWeighted(assignments);
				} else {
					break;
				}
			}
		}
	}

	/**
	 * Calculates a balanced, self-tuning learning rate based on the structural ratio
	 * of target nodes to the square-root compressed color sample pool.
	 *
	 * @param maxNodes The target maximum number of clusters/palette colors (e.g., 32).
	 * @param stdDevSampleSize The total number of elements present in your stdDevSamples list.
	 * @return An optimized learning rate value within a safe stability envelope.
	 */
	public double calculateBalancedLearningRate(int stdDevSampleSize) {
		// Base learning rate anchor for standard Growing Neural Gas transitions
		final double EPSILON_BASE = 0.12; 

		// Safety guard: If the sample list is empty, default safely to the base rate
		if (stdDevSampleSize <= 0) {
			return EPSILON_BASE;
		}

		// 1. Calculate the structural ratio of network capacity to data diversity
		double ratio = (double) maxNodes / stdDevSampleSize;

		// 2. Use square-root scaling to smoothly dampen the adaptation step.
		// This scales down the rate when working with deep color environments, 
		// protecting your nodes from experiencing a catastrophic early hue shift.
		double adaptiveLR = EPSILON_BASE * Math.sqrt(ratio);

		// 3. HARD LIMIT BOUNDARIES (The Stability Guardrails)
		// - Lower bound (0.015) ensures nodes don't freeze up entirely if the pool is massive.
		// - Upper bound (0.080) acts as a hard ceiling. At low target nodes (like 32), 
		//   allowing the learning rate to climb any higher will trigger severe color bleeding.
		return Math.max(0.015, Math.min(0.080, adaptiveLR));
	}

	private GNGNode findBestWinner(double[] sample, List<GNGNode> snapshot) {
		GNGNode winner = null;
		double minDist = Double.MAX_VALUE;
		for (int i = 0; i < snapshot.size(); i++) {
			double d = snapshot.get(i).distance(sample);
			if (d < minDist) {
				minDist = d;
				winner = snapshot.get(i);
			}
		}
		if (winner != null) {
			synchronized(winner) {
				winner.error += minDist;
			}
		}
		return winner;
	}

	public void trainBatch(List<GNGNode> samples, List<GNGNode> uniqueSamples, List<GNGNode> stdDevSamples, int totalEpochs) {
		// 1. PHASE 1: INITIALIZATION
		// Calculate adaptive metrics from your structural equations
		double balancedLR = calculateBalancedLearningRate(stdDevSamples.size());
		int startingPoints = calculateStartingPoints(balancedLR);
		
		// Spread starting nodes across the pure color boundaries
		initializeDistributedNode(uniqueSamples, startingPoints);

		// Split epochs into Growth (Structural) and Tuning (Optimization) stages
		int growthEpochs = (int) (totalEpochs * 0.7); // First 70% of training

		// 2. PHASE 2: TOPOLOGY GROWTH (Using stdDevSamples)
		for (int epoch = 0; epoch < growthEpochs; ++epoch) {
			// Clear historic errors
			for (GNGNode node : nodes)
				node.error = 0.0;
			List<GNGNode> currentNodesSnapshot = new ArrayList<>(nodes);

			// Run matching over the square-root dampened variance stream
			Map<GNGNode, List<GNGNode>> assignments = stdDevSamples.stream()
				.collect(Collectors.groupingBy(sample -> findBestWinner(sample.weight, currentNodesSnapshot)));

			// Adapt node positions using the warm-up smoothing schedule
			double progress = (double) epoch / growthEpochs;
			updateNodeWeightsAdaptive(assignments, balancedLR, progress);

			// Manage ages and insert new nodes until maxNodes budget is met
			manageGraphTopology(assignments, growthEpochs - epoch);
		}

		// 3. PHASE 3: FINAL CENTROID TUNING FOR CRITICAL MSE MINIMIZATION (Using samples)
		// Freeze graph insertions. Force nodes to lock onto true image density.
		int tuningEpochs = totalEpochs - growthEpochs; // Remaining 30% of training
		
		for (int epoch = 0; epoch < tuningEpochs; ++epoch) {
			List<GNGNode> currentNodesSnapshot = new ArrayList<>(nodes);

			// Map the full image pixels directly to the final node layout
			Map<GNGNode, List<GNGNode>> realAssignments = samples.parallelStream()
				.collect(Collectors.groupingByConcurrent(sample -> findBestWinner(sample.weight, currentNodesSnapshot)));

			// STRICT SNAP: Overwrite weights with the true mathematical mean of raw pixels.
			// This step drops the global Mean Squared Error to its absolute floor.
			realAssignments.entrySet().parallelStream().forEach(entry -> {
				GNGNode node = entry.getKey();
				List<GNGNode> cluster = entry.getValue();
				if (cluster == null || cluster.isEmpty())
					return;

				double[] trueMean = new double[node.weight.length];
				for (GNGNode n : cluster) {
					for (int i = 0; i < trueMean.length; ++i) {
						trueMean[i] += n.weight[i];
					}
				}

				// Lock onto the true pixel density center
				for (int i = 0; i < node.weight.length; ++i) {
					node.weight[i] = trueMean[i] / cluster.size();
				}
			});

			// 4. GRAPH TOPOLOGY MANAGEMENT (Ages and Connections)
			// Increment edge age across all current configurations
			for (GNGNode node : nodes) {
				node.incrementAge();
			}

			// Purge old connections that haven't been validated across batches
			for (GNGNode node : nodes) {
				node.removeNeighbourByAge(maxAge);
			}

			// Only remove a node if it is structurally dead AND 
			// failed to capture a single color point during this entire batch epoch.
			nodes.removeIf(node -> {
				if (node.noNeighbor()) {
					List<GNGNode> assigned = realAssignments.get(node);
					return assigned == null || assigned.isEmpty(); 
				}
				return false;
			});

			// 5. DETERMINISTIC GROWTH CALCULATIONS
			int missingNodes = maxNodes - nodes.size();
			if (missingNodes > 0 && epochs - epoch > 0) {
				int targetInsertions = (int) Math.ceil((double) missingNodes / (epochs - epoch));
				for (int i = 0; i < targetInsertions; ++i) {
					if (nodes.size() >= maxNodes)
						break;
					insertNewNodeWeighted(realAssignments);
				}
			}
		}
	}

	private Color[] inxbuild(int nMaxColors) {
		Color[] palette = new Color[nMaxColors];
		int k = 0;
		for (GNGNode node : nodes) {
			double[] channels = node.weight;
			Lab lab1 = new Lab();
			lab1.L = (float) channels[0];
			lab1.A = (float) channels[1];
			lab1.B = (float) channels[2];
			if (channels.length > 3)
				lab1.alpha = (float) channels[3];
			palette[k++] = CIELABConvertor.LAB2RGB(lab1);
			if (k >= nMaxColors)
				break;
		}
		return palette;
	}

	public short nearestColorIndex(final Color[] palette, Color c, final int pos)
	{
		final int offset = getColorIndex(c);
		Short got = nearestMap.get(offset);
		if (got != null)
			return got;

		short k = 0;
		if (c.getAlpha() <= alphaThreshold)
			c = m_transparentColor;
		if(palette.length > 2 && hasAlpha() && c.getAlpha() > alphaThreshold)
			k = 1;
		
		double mindist = 1e100;
		Lab lab1 = getLab(c.getRGB());
		for (short i = k; i < palette.length; ++i) {
			Color c2 = palette[i];

			double curdist = hasSemiTransparency ? BitmapUtilities.sqr(c2.getAlpha() - c.getAlpha()) * TRANS_RATE : 0;
			if (curdist > mindist)
				continue;
			
			Lab lab2 = getLab(c2.getRGB());
			if (palette.length <= 4) {
				curdist = BitmapUtilities.sqr(c2.getRed() - c.getRed())
						+ BitmapUtilities.sqr(c2.getGreen() - c.getGreen()) + BitmapUtilities.sqr(c2.getBlue() - c.getBlue());
				if (hasSemiTransparency)
					curdist += BitmapUtilities.sqr(c2.getAlpha() - c.getAlpha());
			}
			else if (hasSemiTransparency || palette.length < 16) {
				curdist += BitmapUtilities.sqr(lab2.L - lab1.L);
				if (curdist > mindist)
					continue;
				
				curdist += BitmapUtilities.sqr(lab2.A - lab1.A);
				if (curdist > mindist)
					continue;
				
				curdist += BitmapUtilities.sqr(lab2.B - lab1.B);
			}
			else if (palette.length > 32) {
				curdist += Math.abs(lab2.L - lab1.L);
				if (curdist > mindist)
					continue;
				
				curdist += Math.sqrt(BitmapUtilities.sqr(lab2.A - lab1.A) + BitmapUtilities.sqr(lab2.B - lab1.B));
			}
			else {
				double deltaL_prime_div_k_L_S_L = CIELABConvertor.L_prime_div_k_L_S_L(lab1, lab2);
				curdist += BitmapUtilities.sqr(deltaL_prime_div_k_L_S_L);
				if (curdist > mindist)
					continue;
	
				MutableDouble a1Prime = new MutableDouble(), a2Prime = new MutableDouble(), CPrime1 = new MutableDouble(), CPrime2 = new MutableDouble();
				float deltaC_prime_div_k_L_S_L = CIELABConvertor.C_prime_div_k_L_S_L(lab1, lab2, a1Prime, a2Prime, CPrime1, CPrime2);
				curdist += BitmapUtilities.sqr(deltaC_prime_div_k_L_S_L);
				if (curdist > mindist)
					continue;
	
				MutableDouble barCPrime = new MutableDouble(), barhPrime = new MutableDouble();
				float deltaH_prime_div_k_L_S_L = CIELABConvertor.H_prime_div_k_L_S_L(lab1, lab2, a1Prime, a2Prime, CPrime1, CPrime2, barCPrime, barhPrime);
				curdist +=  BitmapUtilities.sqr(deltaH_prime_div_k_L_S_L);
				if (curdist > mindist)
					continue;

				curdist += CIELABConvertor.R_T(barCPrime, barhPrime, deltaC_prime_div_k_L_S_L, deltaH_prime_div_k_L_S_L);
			}

			if (curdist > mindist)
				continue;
			mindist = curdist;
			k = i;
		}
		nearestMap.put(offset, k);
		return k;
	}

	protected short closestColorIndex(final Color[] palette, Color c, final int pos)
	{
		if (c.getAlpha() <= alphaThreshold)
			return nearestColorIndex(palette, c, pos);
		
		final int offset = getColorIndex(c);
		int[] closest = closestMap.get(offset);
		if (closest == null) {
			closest = new int[4];
			closest[2] = closest[3] = Integer.MAX_VALUE;

			for (short k = 0; k < palette.length; ++k) {
				Color c2 = palette[k];
				
				double err = PR * BitmapUtilities.sqr(c2.getRed() - c.getRed());
				if (err >= closest[3])
					continue;

				err += PG * BitmapUtilities.sqr(c2.getGreen() - c.getGreen());
				if (err >= closest[3])
					continue;
				
				err += PB * BitmapUtilities.sqr(c2.getBlue() - c.getBlue());
				if (err >= closest[3])
					continue;

				if(hasSemiTransparency)
					err += PA * BitmapUtilities.sqr(c2.getAlpha() - c.getAlpha());

				if (err < closest[2]) {
					closest[1] = closest[0];
					closest[3] = closest[2];
					closest[0] = k;
					closest[2] = (int) err;
				}
				else if (err < closest[3]) {
					closest[1] = k;
					closest[3] = (int) err;
				}
			}

			if (closest[3] == Integer.MAX_VALUE)
				closest[1] = closest[0];

			closestMap.put(offset, closest);
		}
		
		int idx = 1;
		if (closest[2] == 0 || (random.nextInt(closest[3] + closest[2])) <= closest[3])
			idx = 0;
		
		int MAX_ERR = palette.length;
		if(closest[idx + 2] >= MAX_ERR || closest[idx] == 0 || palette[closest[idx]].getAlpha() < c.getAlpha())
			return nearestColorIndex(palette, c, pos);
		return (short) closest[idx];
	}

	void Clear() {
		m_palette = null;
		saliencies = null;
		closestMap.clear();
		pixelMap.clear();
		nearestMap.clear();
	}

	protected Ditherable getDitherFn(final boolean dither) {
		return new Ditherable() {
			@Override
			public int getColorIndex(Color c) {
				return BitmapUtilities.getColorIndex(c, hasSemiTransparency, hasAlpha());
			}

			@Override
			public short nearestColorIndex(Color[] palette, Color c, final int pos) {
				return DblGNGQuantizer.this.closestColorIndex(palette, c, pos);
			}

		};
	}

	protected short[] dither(Color[] palette, int width, int height, boolean dither)
	{
		Ditherable ditherable = getDitherFn(dither);
		double mDivn = Math.min(0.9, palette.length * 1.0 / pixelMap.size());
		if (hasAlpha())
			mDivn *= -1;

		if (dither && !enforcedDither) {
			short[] qPixels = BitmapUtilities.quantize_image(width, height, pixels, palette, ditherable, hasSemiTransparency, dither);
			Clear();
			return qPixels;
		}

		if (dither) {
			saliencies = new float[pixels.length];
			float saliencyBase = .1f;

			for (int i = 0; i < pixels.length; ++i) {
				Color c = new Color(pixels[i], true);
				Lab lab1 = getLab(c.getRGB());

				saliencies[i] = saliencyBase + (1 - saliencyBase) * lab1.L / 100f * lab1.alpha / 255f;
			}
		}
		
		if (enforcedDither)
			enforcedDither = palette.length < 32 || palette.length > 64;

		short[] qPixels = GilbertCurve.dither(width, height, pixels, palette, ditherable, saliencies, mDivn, dither, enforcedDither);

		if (!dither && palette.length > 32) {
			double delta = mDivn * palette.length;
			mDivn = delta > 0.023 ? 1.0 : (float) (37.013 * delta + 0.906);
			BlueNoise.dither(width, height, pixels, palette, ditherable, qPixels, (float) mDivn);
		}

		Clear();
		return qPixels;
	}

	protected Color[] grabPixels(int[] pixels, int nMaxColors, ThreadLocal<Boolean> isSemiTransparency) {
		final Color[] cPixels = new Color[pixels.length];
		int semiTransCount = 0;
		for (int i = 0; i < pixels.length; ++i) {
			int pixel = pixels[i];
			int alfa = (pixel >> 24) & 0xff;
			cPixels[i] = new Color(pixel, true);
			if (alfa < 0xE0) {
				if (alfa == 0) {
					m_transparentPixelIndex = i;
					if(nMaxColors > 2)
						m_transparentColor = cPixels[i];
					else
						cPixels[i] = m_transparentColor;
				}
				else if(alfa > alphaThreshold)
					++semiTransCount;
			}
		}

		hasSemiTransparency = semiTransCount > 0;
		if(isSemiTransparency != null)
			isSemiTransparency.set(hasSemiTransparency);
		return cPixels;
	}
	
	@FunctionalInterface
	protected interface QuanFn {
		int get(int cnt);
	}
	
	protected Color[] gngquan(final Color[] pixels, int nMaxColors)
	{
		maxNodes = nMaxColors;		// number of colors used
		random = new Random(pixels.length);

		List<GNGNode> samples = new ArrayList<>();
		List<GNGNode> uniqueSamples = new ArrayList<>();
		Map<Integer, Integer> histogram = new HashMap<>();
		for (int i = 0; i < pixels.length; ++i) {
			Color c = BitmapUtilities.getColor(pixels[i], hasSemiTransparency, hasAlpha());
			boolean isRegistered = pixelMap.containsKey(c.getRGB());

			Lab lab1 = getLab(c.getRGB());
			if (hasAlpha())
				samples.add(new GNGNode(new double[]{lab1.L, lab1.A, lab1.B, lab1.alpha}));
			else
				samples.add(new GNGNode(new double[]{lab1.L, lab1.A, lab1.B}));

			if (!isRegistered) {
				if (hasAlpha())
					uniqueSamples.add(new GNGNode(new double[]{lab1.L, lab1.A, lab1.B, lab1.alpha}));
				else
					uniqueSamples.add(new GNGNode(new double[]{lab1.L, lab1.A, lab1.B}));
			}

			Integer got = histogram.get(c.getRGB());
			if (got != null)
				histogram.put(c.getRGB(), got + 1);
			else
				histogram.put(c.getRGB(), 1);
		}

		if(pixelMap.size() <= nMaxColors) {
			/* Fill palette */
			Color[] palette = new Color[pixelMap.size()];
			int k = 0;
			for (Integer pixel : pixelMap.keySet()) {
				Color c = new Color(pixel, hasAlpha());
				palette[k++] = c;

				if(k > 1 && c.getAlpha() == 0) {
					palette[k - 1] = palette[0]; palette[0] = c;
				}
			}

			return palette;
		}

		double mDivn = Math.min(.9, nMaxColors * 1.0 / pixelMap.size());
		List<GNGNode> stdDevSamples = new ArrayList<>();
		for (Integer pixel : histogram.keySet()) {
			int freq = (int) Math.sqrt(histogram.get(pixel));
			for (int i=0; i<freq; ++i) {
				Lab lab1 = getLab(pixel);
				if (hasAlpha())
					stdDevSamples.add(new GNGNode(new double[]{lab1.L, lab1.A, lab1.B, lab1.alpha}));
				else
					stdDevSamples.add(new GNGNode(new double[]{lab1.L, lab1.A, lab1.B}));
			}
		}
		
		if (mDivn < .04 && PG < 1 && PG >= coeffs[0][1] && nMaxColors >= 64)
			enforcedDither = false;
		if (mDivn > .003 && nMaxColors <= 32)
			enforcedDither = false;
		
		if ((nMaxColors < 32 && mDivn > .015 && mDivn < .032) || (nMaxColors >= 32 && nMaxColors < 64 && mDivn > .03 && mDivn < .06))
			trainBatch(uniqueSamples, samples, stdDevSamples, epochs);
		else			
			trainBatch(samples, uniqueSamples, stdDevSamples, epochs);
		
		if (nodes.size() > nMaxColors) {
			System.err.println("Truncated no. of clusters from " + nodes.size() + " to " + nMaxColors);
		}
		else if (nodes.size() < nMaxColors) {
			nMaxColors = nodes.size();
			System.err.println("Reduced no. of clusters to " + nMaxColors);
		}
		return inxbuild(nMaxColors);
	}

	// The work horse for Growing Neural Gas color quantizing.
	public BufferedImage convert(int nMaxColors, boolean dither) {
		if (nMaxColors <= 32)
			PR = PG = PB = PA = 1;
		else {
			PR = coeffs[0][0]; PG = coeffs[0][1]; PB = coeffs[0][2];
		}

		Color[] palette;
		if(m_palette == null) {
			final Color[] cPixels = grabPixels(pixels, nMaxColors, null);
			if (nMaxColors > 2)
				palette = gngquan(cPixels, nMaxColors);
			else {
				palette = new Color[nMaxColors];
				if (hasAlpha()) {
					palette[0] = m_transparentColor;
					palette[1] = Color.BLACK;
				}
				else {
					palette[0] = Color.BLACK;
					palette[1] = Color.WHITE;
				}
			}
		}
		else
			palette = m_palette;

		short[] qPixels = dither(palette, width, height, dither);
		if (hasAlpha() && nMaxColors > 2)
		{
			short k = qPixels[m_transparentPixelIndex];
			palette[k] = m_transparentColor;
		}
		setColorModel(palette);
		return BitmapUtilities.processImagePixels(qPixels, m_colorModel, width, height);
	}

	public Color[] getPalette() {
		return m_palette;
	}

	public boolean hasAlpha() {
		return m_transparentPixelIndex > -1;
	}

	public boolean hasSemiTransparency() {
		return hasSemiTransparency;
	}
}
