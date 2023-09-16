package nQuant.ga;

import java.util.Random;

public interface Chromosome<T extends Chromosome<T> > {

	public float getFitness();

	public T crossover(T mother, int numberOfCrossoverPoints, float crossoverProbability);

	public void mutation(int mutationSize, float mutationProbability);

	public double[] getObjectives();

	public double[] getConvertedObjectives();

	public void resizeConvertedObjectives(int numObj);

	public T makeNewFromPrototype();

	public Random getRandom();

	public boolean dominates(T other);

}
