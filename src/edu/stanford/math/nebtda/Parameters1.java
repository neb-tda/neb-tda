package edu.stanford.math.nebtda;

/**
 * Stores the parameters for 1-cell trials.
 * 
 * @author Henry Adams
 * @author Atanas Atanasov
 */
public class Parameters1 {

	protected double gradConst, convergence, step, nearnessThreshold, clusteringThreshold;
	protected int maxSteps, numTrials, numSuccesses;
	protected InitialBand initialBand;
	protected AnglePenalty anglePenalty;

	/**
	 * Constructor.
	 * 
	 * @param gradConst        		gradient constant
	 * @param convergence       	convergence threshold
	 * @param initialBand			initial band type
	 * @param anglePnFnZero			angle where angle forces first become nonzero
	 * @param anglePnFnOne			angle where angle forces first reach their maximum
	 * @param step		        	size of a step
	 * @param maxSteps				maximum number of steps
	 * @param numTrials				number of trials
	 * @param numSuccesses			number of successes
	 * @param nearnessThreshold		nearness threshold for 0-cells
	 * @param clusteringThreshold	clustering threshold
	 */
	public Parameters1(double gradConst, double convergence, InitialBand initialBand, double anglePnFnZero,
			double anglePnFnOne, double step, int maxSteps, int numTrials, int numSuccesses, double nearnessThreshold,
			double clusteringThreshold) {
		this.gradConst = gradConst;
		this.convergence = convergence;;

		this.initialBand = initialBand;
		this.anglePenalty = new AnglePenaltyTrig(anglePnFnZero, anglePnFnOne);

		this.step = step;
		this.maxSteps = maxSteps;
		this.numTrials = numTrials;
		this.numSuccesses = numSuccesses;

		this.nearnessThreshold = nearnessThreshold;
		this.clusteringThreshold = clusteringThreshold;
	}
	
	/**
	 * Gets the initial band type.
	 * 
	 * @return initial band type
	 */
	public InitialBand getInitialBand() {
		return initialBand;
	}

}
