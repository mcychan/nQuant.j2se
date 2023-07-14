package nQuant.ga;
/*
 * Wu, M.; Yang, D.; Zhou, B.; Yang, Z.; Liu, T.; Li, L.; Wang, Z.; Hu,
 * K. Adaptive Population NSGA-III with Dual Control Strategy for Flexible Job
 * Shop Scheduling Problem with the Consideration of Energy Consumption and Weight. Machines 2021, 9, 344.
 * https://doi.org/10.3390/machines9120344
 * Copyright (c) 2023 Miller Cy Chan
 */

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

public class APNsgaIII<T extends Chromosome<T> > extends NsgaIII<T>
{
	private int _currentGeneration = 0, _max_iterations = 5;
	private int _maxRepeat = Math.min(15, _max_iterations / 2);
	
	// Worst of chromosomes
	protected T _worst;

	// Initializes Adaptive Population NSGA-III with Dual Control Strategy
	public APNsgaIII(T prototype, int numberOfCrossoverPoints, int mutationSize, float crossoverProbability, float mutationProbability)
	{
		super(prototype, numberOfCrossoverPoints, mutationSize, crossoverProbability, mutationProbability);		
	}
	
	private double ex(T chromosome)
	{
		double numerator = 0.0, denominator = 0.0;
		for (int f = 0; f < chromosome.getObjectives().length; ++f) {
			numerator += chromosome.getObjectives()[f] - _best.getObjectives()[f];
			denominator += _worst.getObjectives()[f] - _best.getObjectives()[f];
		}
		return (numerator + 1) / (denominator + 1);
	}
	
	private void popDec(List<T> population)
	{
		int N = population.size();
		if(N <= _populationSize)
			return;
		
		int rank = (int) (.3 * _populationSize);
		
		for(int i = 0; i < N; ++i) {
			double exValue = ex(population.get(i));
			
			if(exValue > .5 && i > rank) {
				population.remove(i);
				if(--N <= _populationSize)
					break;
			}
		}
	}
	
	private void dualCtrlStrategy(List<T> population, int bestNotEnhance, int nMax)
	{
		int N = population.size();
		int nTmp = N;		
		for(int i = 0; i < nTmp; ++i) {
			T chromosome = population.get(i);
			T tumor = chromosome.makeNewFromPrototype();
			tumor.mutation(_mutationSize, _mutationProbability);
			
			_worst = population.get(population.size() - 1);
			if(dominates(tumor, chromosome)) {
				population.set(i, tumor);
				if(dominates(tumor, _best))
					_best = tumor;
			}
			else {
				if(bestNotEnhance >= _maxRepeat && N < nMax) {
					++N;
					if(dominates(_worst, tumor)) {
						population.add(tumor);
						_worst = tumor;
					}
					else
						population.add(population.size() - 1, tumor);
				}
			}
		}
		popDec(population);
	}
	
	@Override
	protected List<T> replacement(List<T> population)
	{
		List<T> result = super.replacement(population);
		result.sort(Comparator.comparing(Chromosome::getFitness, Comparator.reverseOrder()));
		return result;
	}
	
	// Starts and executes algorithm
	@Override
	public void run(int maxRepeat, double minFitness)
    {
		if (_prototype == null)
			return;
		
		List<T>[] pop = new ArrayList[2];
		pop[0] = initialize();
		int nMax = (int) (1.5 * _populationSize);

		int bestNotEnhance = 0;
		double lastBestFit = 0.0;

		int cur = 0, next = 1;
		while(_currentGeneration < _max_iterations)
		{
			T best = getResult();			
			if(_currentGeneration > 0) {	
				double difference = Math.abs(best.getFitness() - lastBestFit);
				if (difference <= 1e-6)
					++bestNotEnhance;
				else {
					lastBestFit = best.getFitness();
					bestNotEnhance = 0;
				}
				
				String status = String.format("\rFitness: %f\t Generation: %d    ", best.getFitness(), _currentGeneration);	
				if(bestNotEnhance >= _maxRepeat)
					status = String.format("\rFitness: %f\t Generation: %d ...", best.getFitness(), _currentGeneration);
				System.out.print(status);
				
				if (best.getFitness() > minFitness) 
					break;

				if (bestNotEnhance > (maxRepeat / 50))		
					reform();
			}				
			
			/******************* crossover *****************/
			List<T> offspring = crossing(pop[cur]);
			
			/******************* mutation *****************/
			offspring.stream().forEach(child -> child.mutation(_mutationSize, _mutationProbability));		

			pop[cur].addAll(offspring);
			
			/******************* replacement *****************/		
			pop[next] = replacement(pop[cur]);			
			_best = dominates(pop[next].get(0), pop[cur].get(0)) ? pop[next].get(0) : pop[cur].get(0);
			
			dualCtrlStrategy(pop[next], bestNotEnhance, nMax);
			
			int temp = cur;
			cur = next;
			next = temp;
			++_currentGeneration;
		}		

	}
	
	@Override
	public String toString()
	{
		return "Adaptive Population NSGA-III with Dual Control Strategy (APNsgaIII)";
	}
}
