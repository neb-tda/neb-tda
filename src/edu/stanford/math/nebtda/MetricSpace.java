package edu.stanford.math.nebtda;

/**
 * Interface for a metric space.
 * 
 * @author Henry Adams
 * @author Atanas Atanasov
 */
public interface MetricSpace<T> {
	
	/**
	 * Computes the distance between two objects.
	 * 
	 * @param a		first object
	 * @param b		second object
	 * @return 		distance
	 */
	double distance(T a, T b);
	
}
