package edu.stanford.math.nebtda;

/**
 * Abstract class for the tangent along a 1-cell.
 * 
 * @author Henry Adams
 * @author Atanas Atanasov
 */
public abstract class Tangent {

	/**
	 * Gets the tangent vector at a node in the band.
	 * 
	 * @param previous		previous node
	 * @param current		current node
	 * @param next			next node
	 * @return				tangent vector at the current node
	 */
	abstract public VectorDouble getTangent(VectorDouble previous, VectorDouble current, VectorDouble next);

}
