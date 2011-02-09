package edu.stanford.math.nebtda;

/**
 * Chooses the representative with the highest value from a cluster.
 * 
 * @author Henry Adams
 * @author Atanas Atanasov
 */
public class ClustRepMaxValue<T> extends ClustRep<T> {

	protected Value<T> valueImplementation;

	/**
	 * Constructor.
	 * 
	 * @param valueImplementation		value class implementation
	 */
	public ClustRepMaxValue(Value<T> valueImplementation) {
		if (valueImplementation == null) {
			throw new IllegalArgumentException("Value implementation reference is null.");
		}
		this.valueImplementation = valueImplementation;
	}

	/**
	 * Gets the representative with the highest value from a cluster.
	 * 
	 * @param cluster		cluster
	 * @return				representative
	 */
	public T getRepresentative(Cluster<T> cluster) {
		
		if (cluster.isEmpty()) {
			throw new IllegalArgumentException("The cluster is empty.");
		}
		T max = cluster.removeFirst();
		T cur;
		while(!cluster.isEmpty()) {
			cur = cluster.removeFirst();
			if (valueImplementation.getValue(cur) > valueImplementation.getValue(max)) {
				max = cur;
			}
		}
		return max;
	}

}
