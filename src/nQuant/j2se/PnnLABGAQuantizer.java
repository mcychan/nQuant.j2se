package nQuant.j2se;

/* Fast pairwise nearest neighbor based genetic algorithm with CIELAB color space genetic algorithm
* Copyright (c) 2023 Miller Cy Chan */

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.util.Arrays;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.ConcurrentHashMap;

import nQuant.ga.Chromosome;
import nQuant.j2se.CIELABConvertor.Lab;

public class PnnLABGAQuantizer implements AutoCloseable, Chromosome<PnnLABGAQuantizer> {
	private double _fitness = Double.NEGATIVE_INFINITY;
	private double ratioX = 0, ratioY = 0;
	private double[] _convertedObjectives;
	private double[] _objectives;
	private Color[] cPixels;
	private Random _random;
	private PnnLABQuantizer m_pq;
	
	private static int _dp = 1, _nMaxColors = 256;
	private static double minRatio = 0, maxRatio = 1.0;

	private static Map<String, double[]> _fitnessMap = new ConcurrentHashMap<>();

	public PnnLABGAQuantizer(PnnLABQuantizer pq, int nMaxColors) {
		// increment value when criteria violation occurs
		_objectives = new double[4];
		_random = new Random(pq.getPixels().length);
		m_pq = new PnnLABQuantizer(pq);
		if(pq.isGA()) {
			cPixels = m_pq.getCPixels();
			return;
		}		
		
		_nMaxColors = nMaxColors;
		
		int[] pixels = m_pq.getPixels();
		ThreadLocal<Boolean> isSemiTransparency = new ThreadLocal<>();
		cPixels = m_pq.grabPixels(pixels, nMaxColors, isSemiTransparency);
		boolean hasSemiTransparency = isSemiTransparency.get();
		minRatio = (hasSemiTransparency || nMaxColors < 64) ? .0111 : .85;
		maxRatio = Math.min(1.0, nMaxColors / ((nMaxColors < 64) ? 500.0 : 50.0));
		_dp = maxRatio < .1 ? 10000 : 100;
	}	

	public double[] getRatios() {
		return new double[] {ratioX, ratioY};
	}
	
	private String getRatioKey() {
		StringBuffer sb = new StringBuffer();
		sb.append((int) (ratioX * _dp));
		double difference = Math.abs(ratioX - ratioY);
		if (difference <= 0.0000001)
			return sb.toString();
		
		sb.append(";").append((int) (ratioY * _dp * 100));
		return sb.toString();
	}
	
	@FunctionalInterface
	protected interface AmplifyFn {
		double get(double val);
	}
	
	private AmplifyFn getAmplifyFn(boolean tooSmall) {
		if (tooSmall)
			return val -> Math.exp(val);
		return val -> Math.log(val);
	}
	
	private void calculateError(double[] errors) {
		double maxError = maxRatio < .1 ? .5 : .0625;
		if(m_pq.hasAlpha())
			maxError = 1;
		
		double fitness = 0;
		boolean tooSmall = false; // any error < exp(1.0) concluded as too small
		for(int i = 0; i < errors.length; ++i) {
			errors[i] /= maxError * cPixels.length;
			if(errors[i] < 3)
				tooSmall = true;
		}		
		
		AmplifyFn amplifyFn = getAmplifyFn(tooSmall);
		for(int i = 0; i < errors.length; ++i) {
			if(i == 0 && errors[i] > maxError)
				errors[i] *= amplifyFn.get(errors[i]);
			else if(errors[i] > 2 * maxError)
				errors[i] *= amplifyFn.get(errors[i]);
			fitness -= errors[i];
		}
		
		_objectives = errors;
		_fitness = fitness;
	}
	
	private void calculateFitness() {
		final String ratioKey = getRatioKey();
		if(_fitnessMap.containsKey(ratioKey)) {
			_objectives = _fitnessMap.get(ratioKey);
			_fitness = -1 * Arrays.stream(_objectives).sum();
			return;
		}	
		
		m_pq.setRatio(ratioX, ratioY);
		Color[] palette = m_pq.pnnquan(cPixels, _nMaxColors);

		int threshold = maxRatio < .1 ? -64 : -112;
		double[] errors = new double[_objectives.length];
		for (int i = 0; i < cPixels.length; ++i) {
			if(BlueNoise.RAW_BLUE_NOISE[i & 4095] > threshold)
				continue;
			
			Color c = cPixels[i];
			Lab lab1 = m_pq.getLab(cPixels[i].getRGB());
			short qPixelIndex = m_pq.nearestColorIndex(palette, c, i);
			Lab lab2 = m_pq.getLab(palette[qPixelIndex].getRGB());
						
			if(m_pq.hasAlpha()) {
				errors[0] += BitmapUtilities.sqr(lab2.L - lab1.L);
				errors[1] += BitmapUtilities.sqr(lab2.A - lab1.A);
				errors[2] += BitmapUtilities.sqr(lab2.B - lab1.B);
				errors[3] += BitmapUtilities.sqr(lab2.alpha - lab1.alpha) / Math.exp(1.5);
			}
			else {
				errors[0] += Math.abs(lab2.L - lab1.L);
				errors[1] += Math.sqrt(BitmapUtilities.sqr(lab2.A - lab1.A) + BitmapUtilities.sqr(lab2.B - lab1.B));
			}
		}
		calculateError(errors);
		_fitnessMap.put(ratioKey, _objectives);
	}

	@Override
	public void close() throws Exception {
		cPixels = null;
		_fitnessMap.clear();
	}
	
	private double randrange(double min, double max)
	{
		return min + _random.nextDouble() * (max - min);
	}

	public void setRatio(double ratioX, double ratioY) {
		double difference = Math.abs(ratioX - ratioY);
		if (difference <= minRatio)
			ratioY = ratioX;
		this.ratioX = Math.min(Math.max(ratioX, minRatio), maxRatio);
		this.ratioY = Math.min(Math.max(ratioY, minRatio), maxRatio);
	}

	@Override
	public float getFitness() {
		return (float) _fitness;
	}
	
	private double rotateLeft(double u, double v, double delta) {
		double theta = Math.PI * randrange(minRatio, maxRatio) / Math.exp(delta);
		double result = u * Math.sin(theta) + v * Math.cos(theta);
		if(result <= minRatio || result >= maxRatio)
			result = rotateLeft(u, v, delta + .5);
		return result;
	}
	
	private double rotateRight(double u, double v, double delta) {
		double theta = Math.PI * randrange(minRatio, maxRatio) / Math.exp(delta);
		double result = u * Math.cos(theta) - v * Math.sin(theta);
		if(result <= minRatio || result >= maxRatio)
			result = rotateRight(u, v, delta + .5);
		return result;
	}

	@Override
	public PnnLABGAQuantizer crossover(PnnLABGAQuantizer mother, int numberOfCrossoverPoints, float crossoverProbability) {		
		PnnLABGAQuantizer child = makeNewFromPrototype();
		if (_random.nextInt(100) <= crossoverProbability)
			return child;
		
		double ratioX = rotateRight(this.ratioX, mother.getRatios()[1], 0.0);
		double ratioY = rotateLeft(this.ratioY, mother.getRatios()[0], 0.0);
		child.setRatio(ratioX, ratioY);
		child.calculateFitness();
		return child;
	}
	
	private double boxMuller(double value) {
		double r1 = randrange(minRatio, maxRatio);
		return Math.sqrt(-2 * Math.log(value)) * Math.cos(2 * Math.PI * r1);
	}

	@Override
	public void mutation(int mutationSize, float mutationProbability) {
		// check probability of mutation operation
		if (_random.nextInt(100) > mutationProbability)
			return;
		
		double ratioX = this.ratioX;
		double ratioY = this.ratioY;
		if(_random.nextDouble() > .5)
			ratioX = boxMuller(ratioX);
		else
			ratioY = boxMuller(ratioY);
		
		setRatio(ratioX, ratioY);
		calculateFitness();
	}

	@Override
	public double[] getObjectives() {
		return _objectives;
	}

	@Override
	public double[] getConvertedObjectives() {
		return _convertedObjectives;
	}

	@Override
	public void resizeConvertedObjectives(int numObj) {
		_convertedObjectives = new double[numObj];
	}

	@Override
	public PnnLABGAQuantizer makeNewFromPrototype() {
		PnnLABGAQuantizer child = new PnnLABGAQuantizer(m_pq, _nMaxColors);
		double minRatio2 = 2 * minRatio;
		if(minRatio2 > 1)
			minRatio2 = 0;
		double ratioX = randrange(minRatio, maxRatio);
		double ratioY = ratioX < minRatio2 ? randrange(minRatio, maxRatio) : ratioX;		
		child.setRatio(ratioX, ratioY);
		child.calculateFitness();
		return child;
	}
	
	public Random getRandom() {
		return _random;
	}
	
	public BufferedImage convert(boolean dither) {
		m_pq.setRatio(ratioX, ratioY);
		m_pq.pnnquan(cPixels, _nMaxColors);
		return m_pq.convert(_nMaxColors, dither);
	}

	public boolean hasAlpha() {
		return m_pq.hasAlpha();
	}

	public String getResult() {
		double difference = Math.abs(ratioX - ratioY);
		if (difference <= 0.0000001)
			return String.valueOf(ratioX);
		return ratioX + ", " + ratioY;
	}

}
