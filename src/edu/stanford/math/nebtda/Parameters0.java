package edu.stanford.math.nebtda;

/**
 * Stores the parameters for 0-cell trials.
 * 
 * @author Henry Adams
 * @author Atanas Atanasov
 */
public class Parameters0 {
	
	protected double convergence;
	protected int numTrials;
	protected double clusteringThreshold;	

	/**
	 * Constructor.
	 * 
	 * @param convergence				convergence threshold
	 * @param numTrials					number of trials
	 * @param clusteringThreshold		clustering threshold
	 */
	public Parameters0(double convergence, int numTrials, double clusteringThreshold) {
		this.convergence = convergence;
		this.numTrials = numTrials;
		this.clusteringThreshold = clusteringThreshold;
	}

}
