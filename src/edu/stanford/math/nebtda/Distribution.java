package edu.stanford.math.nebtda;

/**
 * Abstract class which declares the functionality of a distribution.
 * 
 * @author Henry Adams
 * @author Atanas Atanasov
 */
public abstract class Distribution {
	
	/**
	 * Empty constructor.
	 */
	public Distribution() {
	}

	/**
	 * Computes the value of the distribution at a point.
	 * 
	 * @param point		point
	 * @param mean		mean of the distribution
	 * @param std		standard deviation of the distribution
	 * @return			value at the point
	 */
	public abstract double getValue(VectorDouble point, VectorDouble mean, double std);

	/**
	 * Computes the gradient of the distribution at a point.
	 * 
	 * @param point		point
	 * @param mean		mean of the distribution
	 * @param std		standard deviation of the distribution
	 * @return			gradient at the point
	 */
	public abstract VectorDouble getGradient(VectorDouble point, VectorDouble mean, double std);

	/**
	 * Computes the value and gradient of the distribution at a point.
	 * 
	 * @param point		point
	 * @param mean		mean of the distribution
	 * @param std		standard deviation of the distribution
	 * @return			value and gradient at the point
	 */
	public abstract Tuple<Double, VectorDouble> getValueGradient(VectorDouble point, VectorDouble mean, double std);
	
	/**
	 * Computes the maximum norm of the gradient of the distribution.
	 * 
	 * @param std			standard deviation of the distribution
	 * @param dimension		dimension of the distribution
	 * @return				maximum norm of the gradient
	 */
	public abstract double maxGradNorm(double std, int dimension);
}