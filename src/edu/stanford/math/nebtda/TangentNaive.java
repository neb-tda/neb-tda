package edu.stanford.math.nebtda;

/**
 * Implements the naive tangent along a 1-cell, which is given by averaging the two adjacent edges.
 * 
 * @author Henry Adams
 * @author Atanas Atanasov
 */
public class TangentNaive extends Tangent {

	/**
	 * Empty constructor.
	 */
	public TangentNaive() {
	}

	/**
	 * Gets the tangent vector at a node in the band.
	 * 
	 * @param previous		previous node
	 * @param current		current node
	 * @param next			next node
	 * @return				tangent vector at the current node
	 */
	public VectorDouble getTangent(VectorDouble previous, VectorDouble current, VectorDouble next) {
		return VectorDouble.subtract(next, previous);
	}

}
