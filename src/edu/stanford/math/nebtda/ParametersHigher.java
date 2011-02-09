package edu.stanford.math.nebtda;

/**
 * Stores the parameters for a higher dimensional cell trial.
 * 
 * @author Henry Adams
 * @author Atanas Atanasov
 */
public class ParametersHigher {
	
	protected double gradConst, convergence, step;
	protected int maxSteps;
	
	/**
	 * Constructor.
	 * 
	 * @param gradConst        		gradient constant
	 * @param convergence       	convergence threshold
	 * @param step		        	size of a step
	 * @param maxSteps				maximum number of steps
	 */
	public ParametersHigher(double gradConst, double convergence, double step, int maxSteps) {
		this.gradConst = gradConst;
		this.convergence = convergence;
		this.step = step;
		this.maxSteps = maxSteps;
	}

}
