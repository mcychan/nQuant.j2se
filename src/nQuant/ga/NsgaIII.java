package nQuant.ga;

/*
 * Deb K , Jain H . An Evolutionary Many-Objective Optimization Algorithm Using Reference Point-Based Nondominated Sorting Approach,
 * Part I: Solving Problems With Box Constraints[J]. IEEE Transactions on Evolutionary Computation, 2014, 18(4):577-601.
 * Copyright (c) 2023 Miller Cy Chan
 */

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

public class NsgaIII<T extends Chromosome<T> >
{
	// Best of chromosomes
	protected T _best;

	// Prototype of chromosomes in population
	protected T _prototype;
	
	// Number of chromosomes
	protected int _populationSize;

	// Number of crossover points of parent's class tables
	protected int _numberOfCrossoverPoints;

	// Number of classes that is moved randomly by single mutation operation
	protected int _mutationSize;

	// Probability that crossover will occur
	protected float _crossoverProbability;

	// Probability that mutation will occur
	protected float _mutationProbability;
	
	protected List<Integer> _objDivision;
	
	private int _criteriaLength;
	
	private static Random _random;

	// Initializes NsgaIII
	private NsgaIII(T prototype, int numberOfChromosomes)
	{
		_prototype = prototype;
		_random = prototype.getRandom();
		
		// there should be at least 2 chromosomes in population
		if (numberOfChromosomes < 2)
			numberOfChromosomes = 2;
		_populationSize = numberOfChromosomes;
	}

	public NsgaIII(T prototype, int numberOfCrossoverPoints, int mutationSize, float crossoverProbability, float mutationProbability)
	{
		this(prototype, 9);
		_criteriaLength = prototype.getObjectives().length;
		_mutationSize = mutationSize;
		_numberOfCrossoverPoints = numberOfCrossoverPoints;
		_crossoverProbability = crossoverProbability;
		_mutationProbability = mutationProbability;
		
		_objDivision = new ArrayList<>();
		if(_criteriaLength < 8)
			_objDivision.add(6);
		else {
			_objDivision.add(3);
			_objDivision.add(2);
		}
	}

	// Returns pointer to best chromosomes in population
	public T getResult()
	{
		return _best;
	}
	
	private static int rand(int size)
	{
		return _random.nextInt(size);
	}
	
	private static class ReferencePoint {
		private int memberSize;
		private double[] position;
		private Map<Integer, Double> potentialMembers;
		
		ReferencePoint(int M) {
			memberSize = 0;
			position = new double[M];
			potentialMembers = new HashMap<>();
		}
		
		static void generateRecursive(List<ReferencePoint> rps, ReferencePoint pt, int numObjs, int left, int total, int element) {
			if (element == numObjs - 1) {
				pt.position[element] = left * 1.0 / total;
				rps.add(pt);
			}
			else {
				for (int i = 0; i <= left; ++i) {
					pt.position[element] = i * 1.0 / total;
					generateRecursive(rps, pt, numObjs, left - i, total, element + 1);
				}
			}
		}
		
		void addMember()
		{
			++memberSize;
		}
		
		void addPotentialMember(Integer memberInd, double distance)
		{
			Double currDistance = potentialMembers.get(memberInd);
			if(currDistance == null || distance < currDistance)
				potentialMembers.put(memberInd, distance);
		}
		
		int findClosestMember()
		{
			return Collections.min(potentialMembers.entrySet(), Map.Entry.comparingByValue()).getKey();
		}
		
		boolean hasPotentialMember()
		{
			return !potentialMembers.isEmpty();
		}

		int randomMember()
		{
			if (potentialMembers.isEmpty())
				return -1;

			Integer[] members = potentialMembers.keySet().toArray(new Integer[0]);
			return members[rand(potentialMembers.size())];			
		}
		
		void removePotentialMember(Integer memberInd)
		{
			potentialMembers.remove(memberInd);
		}

		static void generateReferencePoints(List<ReferencePoint> rps, int M, final List<Integer> p) {
			ReferencePoint pt = new ReferencePoint(M);
			generateRecursive(rps, pt, M, p.get(0), p.get(0), 0);

			if (p.size() > 1) { // two layers of reference points (Check Fig. 4 in NSGA-III paper)
				List<ReferencePoint> insideRps = new ArrayList<>();
				generateRecursive(insideRps, pt, M, p.get(1), p.get(1), 0);

				double center = 1.0 / M;

				for (ReferencePoint insideRp : insideRps) {
					for (int j = 0; j < insideRp.position.length; ++j)
						insideRp.position[j] = center + insideRp.position[j] / 2; // (k=num_divisions/M, k, k, ..., k) is the center point

					rps.add(insideRp);
				}
			}
		}
		
	}
	
	private static double perpendicularDistance(final double[] direction, final double[] point)
	{
		double numerator = 0, denominator = 0;
		for (int i = 0; i < direction.length; ++i) {
			numerator += direction[i] * point[i];
			denominator += Math.pow(direction[i], 2);
		}
		
		if(denominator <= 0)
			return Double.MAX_VALUE;
		
		double k = numerator / denominator;
		double d = 0;
		for (int i = 0; i < direction.length; ++i)
			d += Math.pow(k * direction[i] - point[i], 2);

		return Math.sqrt(d);
	}
	
	private void associate(List<ReferencePoint> rps, final List<T> pop, final List<List<Integer> > fronts) {
		for (int t = 0; t < fronts.size(); ++t) {
			for (Integer memberInd : fronts.get(t)) {
				int minRp = rps.size() - 1;
				double minDist = Double.MAX_VALUE;
				for (int r = 0; r < rps.size(); ++r) {
					double d = perpendicularDistance(rps.get(r).position, pop.get(memberInd).getConvertedObjectives());
					if (d < minDist) {
						minDist = d;
						minRp = r;
					}
				}

				if (t + 1 != fronts.size()) // associating members in St/Fl (only counting)
					rps.get(minRp).addMember();
				else
					rps.get(minRp).addPotentialMember(memberInd, minDist);

			}// for - members in front
		}// for - fronts
	}
	

	private static double[] guassianElimination(List<Double>[] A, final double[] b)
	{
		final int N = A.length;
		for (int i = 0; i < N; ++i)
			A[i].add(b[i]);

		for (int base = 0; base < N - 1; ++base) {
			for (int target = base + 1; target < N; ++target) {
				double ratio = A[target].get(base) / A[base].get(base);
				for (int term = 0; term < A[base].size(); ++term)
					A[target].set(term, A[target].get(term) - A[base].get(term) * ratio);
			}
		}

		double[] x = new double[N];
		for (int i = N - 1; i >= 0; --i) {
			for (int known = i + 1; known < N; ++known)
				A[i].set(N, A[i].get(N) - A[i].get(known) * x[known]);

			x[i] = A[i].get(N) / A[i].get(i);
		}
		return x;
	}
	
	// ----------------------------------------------------------------------
	// ASF: Achivement Scalarization Function
	// ----------------------------------------------------------------------
	private static double ASF(final double[] objs, final double[] weight)
	{
		double max_ratio = -Double.MAX_VALUE;
		for (int f = 0; f < objs.length; ++f) {
			double w = Math.max(weight[f], 1e-6);
			max_ratio = Math.max(max_ratio, objs[f] / w);
		}
		return max_ratio;
	}
	

	private List<Integer> findExtremePoints(final List<T> pop, final List<List<Integer> > fronts) {
		final int numObj = pop.get(0).getObjectives().length;
		
		List<Integer> exp = new ArrayList<>();
		for (int f = 0; f < numObj; ++f) {
			double[] w = DoubleStream.generate(() -> 1e-6).limit(numObj).toArray();
			w[f] = 1.0;

			double minASF = Double.MAX_VALUE;
			int minIndv = fronts.get(0).size();

			for (int frontIndv : fronts.get(0)) { // only consider the individuals in the first front
				double asf = ASF(pop.get(frontIndv).getConvertedObjectives(), w);

				if (asf < minASF) {
					minASF = asf;
					minIndv = frontIndv;
				}
			}

			exp.add(minIndv);
		}

		return exp;
	}
	

	private double[] findMaxObjectives(final List<T> pop)
	{
		final int numObj = pop.get(0).getObjectives().length;
		double[] maxPoint = DoubleStream.generate(() -> -Double.MAX_VALUE).limit(numObj).toArray();
		for (int i = 0; i < pop.size(); ++i) {
			for (int f = 0; f < maxPoint.length; ++f)
				maxPoint[f] = Math.max(maxPoint[f], pop.get(i).getObjectives()[f]);
		}

		return maxPoint;
	}
	

	private static int findNicheReferencePoint(final List<ReferencePoint> rps)
	{
		// find the minimal cluster size
		int minSize = Integer.MAX_VALUE;
		for (ReferencePoint rp : rps)
			minSize = Math.min(minSize, rp.memberSize);

		// find the reference points with the minimal cluster size Jmin
		List<Integer> minRps = new ArrayList<>();
		for (int r = 0; r < rps.size(); ++r) {
			if (rps.get(r).memberSize == minSize)
				minRps.add(r);
		}

		// return a random reference point (j-bar)
		return minRps.get(rand(minRps.size()));
	}
	

	private List<Double> constructHyperplane(final List<T> pop, final List<Integer> extremePoints)
	{
		final int numObj = pop.get(0).getObjectives().length;
		// Check whether there are duplicate extreme points.
		// This might happen but the original paper does not mention how to deal with it.
		boolean duplicate = false;
		for (int i = 0; !duplicate && i < extremePoints.size(); ++i) {
			for (int j = i + 1; !duplicate && j < extremePoints.size(); ++j)
				duplicate = (extremePoints.get(i) == extremePoints.get(j));
		}

		List<Double> intercepts = new ArrayList<>();

		boolean negativeIntercept = false;
		if (!duplicate) {
			// Find the equation of the hyperplane
			double[] b = DoubleStream.generate(() -> 1.0).limit(numObj).toArray();
			List<Double>[] A = new ArrayList[extremePoints.size()];
			for (int p = 0; p < extremePoints.size(); ++p)
			{
				double[] convertedObjs = pop.get( extremePoints.get(p) ).getConvertedObjectives();
				A[p] = DoubleStream.of(convertedObjs).boxed().collect(Collectors.toList());
			}
			
			double[] x = guassianElimination(A, b);
			// Find intercepts
			for (int f = 0; f < numObj; ++f) {
				intercepts.add(1.0 / x[f]);

				if(x[f] < 0) {
					negativeIntercept = true;
					break;
				}
			}
		}

		if (duplicate || negativeIntercept) { // follow the method in Yuan et al. (GECCO 2015)
			double[] maxObjs = findMaxObjectives(pop);
			intercepts = DoubleStream.of(maxObjs).boxed().collect(Collectors.toList());
		}
		return intercepts;
	}
	

	private void normalizeObjectives(List<T> pop, final List<List<Integer> > fronts, final List<Double> intercepts, final List<Double> idealPoint)
	{		
		for (List<Integer> front : fronts) {
			for (int ind : front) {
				double[] convObjs = pop.get(ind).getConvertedObjectives();
				for (int f = 0; f < convObjs.length; ++f) {
					if (Math.abs(intercepts.get(f) - idealPoint.get(f)) > 10e-10) // avoid the divide-by-zero error
						convObjs[f] /= intercepts.get(f) - idealPoint.get(f);
					else
						convObjs[f] /= 10e-10;
				}
			}
		}

	}
	
	protected List<List<Integer> > nondominatedSort(List<T> pop) {
		List<List<Integer> > fronts = new ArrayList<>();
		int numAssignedIndividuals = 0;
		int rank = 1;
		int[] indvRanks = new int[pop.size()];

		while (numAssignedIndividuals < pop.size()) {
			List<Integer> curFront = new ArrayList<>();

			for (int i = 0; i < pop.size(); ++i) {
				if (indvRanks[i] > 0)
					continue; // already assigned a rank

				boolean beDominated = false;
				for (int j = 0; j < curFront.size(); ++j) {
					if (pop.get( curFront.get(j) ).dominates(pop.get(i)) ) { // i is dominated
						beDominated = true;
						break;
					}
					else if (pop.get(i).dominates( pop.get( curFront.get(j) )) ) // i dominates a member in the current front
						curFront.remove(j--);
				}
				
				if (!beDominated)
					curFront.add(i);
			}

			for (int front : curFront)
				indvRanks[front] = rank;

			fronts.add(curFront);
			numAssignedIndividuals += curFront.size();
			
			++rank;
		}

		return fronts;
	}
	

	private static int selectClusterMember(final ReferencePoint rp) {
		if (rp.hasPotentialMember()) {
			if (rp.memberSize == 0) // currently has no member
				return rp.findClosestMember();

			return rp.randomMember();
		}
		return -1;
	}
	

	private List<Double> translateObjectives(List<T> pop, final List<List<Integer> > fronts)
	{
		List<Double> idealPoint = new ArrayList<>();
		final int numObj = pop.get(0).getObjectives().length;
		for (int f = 0; f < numObj; ++f) {
			double minf = Double.MAX_VALUE;
			for (int frontIndv : fronts.get(0)) // min values must appear in the first front
				minf = Math.min(minf, pop.get(frontIndv).getObjectives()[f]);
			
			idealPoint.add(minf);

			for (List<Integer> front : fronts) {
				for (int ind : front) {
					T chromosome = pop.get(ind);
					chromosome.resizeConvertedObjectives(numObj);
					chromosome.getConvertedObjectives()[f] = chromosome.getObjectives()[f] - minf;
				}
			}
		}
		
		return idealPoint;
	}

	private List<T> selection(List<T> cur, List<ReferencePoint> rps) {
		List<T> next = new ArrayList<>();
		
		// ---------- Step 4 in Algorithm 1: non-dominated sorting ----------
		List<List<Integer> > fronts = nondominatedSort(cur);
		
		// ---------- Steps 5-7 in Algorithm 1 ----------
		int last = 0, next_size = 0;
		while (next_size < _populationSize) {
			next_size += fronts.get(last++).size();
		}
		
		fronts = fronts.subList(0, last); // remove useless individuals

		for (int t = 0; t < fronts.size() - 1; ++t) {
			for (int frontIndv : fronts.get(t))
				next.add(cur.get(frontIndv));
		}

		// ---------- Steps 9-10 in Algorithm 1 ----------
		if (next.size() == _populationSize)
			return next;

		// ---------- Step 14 / Algorithm 2 ----------
		List<Double> idealPoint = translateObjectives(cur, fronts);
		
		List<Integer> extremePoints = findExtremePoints(cur, fronts);

		List<Double> intercepts = constructHyperplane(cur, extremePoints);

		normalizeObjectives(cur, fronts, intercepts, idealPoint);

		// ---------- Step 15 / Algorithm 3, Step 16 ----------
		associate(rps, cur, fronts);

		// ---------- Step 17 / Algorithm 4 ----------
		while (next.size() < _populationSize) {
			int minRp = findNicheReferencePoint(rps);

			int chosen = selectClusterMember(rps.get(minRp));
			if (chosen < 0) // no potential member in Fl, disregard this reference point
				rps.remove(minRp);
			else {
				rps.get(minRp).addMember();
				rps.get(minRp).removePotentialMember(chosen);
				next.add(cur.get(chosen));
			}
		}

		return next;
	}

	protected List<T> crossing(List<T> population)
	{
		int populationSize = population.size();
		List<T> offspring = new ArrayList<T>(populationSize);
		IntStream.range(0, populationSize).parallel().forEach(i -> {
			if(i % 2 == 0) {
				int father = rand(populationSize), mother = rand(populationSize);				
				offspring.add(population.get(father).crossover(population.get(mother), _numberOfCrossoverPoints, _crossoverProbability));
				if((i + 1) < populationSize)
					offspring.add(population.get(mother).crossover(population.get(father), _numberOfCrossoverPoints, _crossoverProbability));
			}
		});
		return offspring;
	}

	protected List<T> initialize()
	{
		List<T> result = new ArrayList<T>(_populationSize);
		result.add(_prototype.makeNewFromPrototype());
		// initialize new population with chromosomes randomly built using prototype
		IntStream.range(1, _populationSize).parallel()
			.forEach(i -> result.add(_prototype.makeNewFromPrototype()));
		return result;
	}
	
	protected void reform()
	{
		_random = new Random(System.currentTimeMillis());
		if(_crossoverProbability < 95)
			_crossoverProbability += 1.0f;
		else if(_mutationProbability < 30)
			_mutationProbability += 1.0f;
	}
	
	protected List<T> replacement(List<T> population)
	{
		List<ReferencePoint> rps = new ArrayList<>();
		ReferencePoint.generateReferencePoints(rps, _criteriaLength, _objDivision);			
		return selection(population, rps);
	}
	
	// Starts and executes algorithm
	public void run(int maxRepeat, double minFitness)
	{
		if (_prototype == null)
			return;		
		
		List<T>[] pop = new ArrayList[2];
		pop[0] = initialize();

		// Current generation
		int currentGeneration = 0;
		int bestNotEnhance = 0;
		double lastBestFit = 0.0;

		int cur = 0, next = 1;
		for (; ;)
		{				
			T best = getResult();
			if(currentGeneration > 0) {
				String status = String.format("\rFitness: %f\t Generation: %d", best.getFitness(), currentGeneration);
				System.out.print(status);
				
				// algorithm has reached criteria?
				if (best.getFitness() > minFitness)
					break;
	
				double difference = Math.abs(best.getFitness() - lastBestFit);
				if (difference <= 0.0000001)
					++bestNotEnhance;
				else {
					lastBestFit = best.getFitness();
					bestNotEnhance = 0;
				}

				if (bestNotEnhance > (maxRepeat / 50))
					reform();
			}
			
			/******************* crossover *****************/
			List<T> offspring = crossing(pop[cur]);
			
			/******************* mutation *****************/
			offspring.parallelStream().forEach(child -> child.mutation(_mutationSize, _mutationProbability));
			
			pop[cur].addAll(offspring);
			
			/******************* replacement *****************/	
			pop[next] = replacement(pop[cur]);
			_best = pop[next].get(0).dominates( pop[cur].get(0)) ? pop[next].get(0) : pop[cur].get(0);
			
			int temp = cur;
			cur = next;
			next = temp;
			++currentGeneration;
		}
	}	

	@Override
	public String toString()
	{
		return "NSGA III";
	}
}

