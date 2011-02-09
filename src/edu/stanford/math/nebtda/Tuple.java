package edu.stanford.math.nebtda;

/**
 * Implements an ordered pair.
 * 
 * @author Henry Adams
 * @author Atanas Atanasov
 */
public final class Tuple<T1, T2> {
	
	protected T1 obj1;
	protected T2 obj2;

	/**
	 * Constructor.
	 * 
	 * @param obj1		first object
	 * @param obj2		second object
	 */
	public Tuple(T1 obj1, T2 obj2) {
		this.obj1 = obj1;
		this.obj2 = obj2;
	}

	/**
	 * Gets the first object.
	 * 
	 * @return first object
	 */
	public T1 get1() {
		return obj1;
	}

	/**
	 * Gets the second object.
	 * 
	 * @return second object
	 */
	public T2 get2() {
		return obj2;
	}
}
