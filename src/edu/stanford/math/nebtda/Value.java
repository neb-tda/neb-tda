package edu.stanford.math.nebtda;

/**
 * Interface for objects which can return a value, typically density.
 * 
 * @author Henry Adams
 * @author Atanas Atanasov
 */
public interface Value<T> {
	
	/**
	 * Gets the value.
	 * 
	 * @param t		object
	 * @return		value of the object
	 */
	double getValue(T t);
	
}