package edu.stanford.math.nebtda;

/**
 * Implements the zero angle penalty on a 1-cell.
 * 
 * @author Henry Adams
 * @author Atanas Atanasov
 */
public class AnglePenaltyNone extends AnglePenalty {

	/**
	 * Empty constructor.
	 */
	public AnglePenaltyNone() {
	}

	/**
	 * Computes the angle penalty at a node, which is zero.
	 * 
	 * @param u_plus		first adjacent edge
	 * @param u_minus		second adjacent edge
	 * @return 				the zero vector
	 */
	public VectorDouble getAnglePenalty(VectorDouble u_plus, VectorDouble u_minus) {
		return new VectorDouble(u_plus.coordinates.length);
	}
	
}
