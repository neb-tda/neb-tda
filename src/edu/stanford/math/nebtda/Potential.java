package edu.stanford.math.nebtda;

/**
 * Abstract class for a differentiable potential function.
 * 
 * @author Henry Adams
 * @author Atanas Atanasov
 */
abstract public class Potential {
	
	/**
	 * Computes the value at a point.
	 * 
	 * @param point		point
	 * @return 			value
	 */
	public abstract double getValue(VectorDouble point);
	
	/**
	 * Computes the gradient vector at a point.
	 * 
	 * @param point		point
	 * @return 			gradient vector
	 */
	public abstract VectorDouble getGradient(VectorDouble point);

	/**
	 * Computes the value and gradient vector at a point.
	 * 
	 * @param point		point
	 * @return 			value and gradient vector
	 */
	public abstract Tuple<Double, VectorDouble> getValueGradient(VectorDouble point);

}