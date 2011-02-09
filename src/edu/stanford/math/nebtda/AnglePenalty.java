package edu.stanford.math.nebtda;

/**
 * Abstract class for the angle penalty on a 1-cell, which acts as a smoothing force.
 * 
 * @author Henry Adams
 * @author Atanas Atanasov
 */
public abstract class AnglePenalty {

	/**
	 * Computes the angle penalty at a node.
	 * 
	 * @param u_plus		first adjacent edge
	 * @param u_minus		second adjacent edge
	 * @return 				angle penalty
	 */
	abstract public VectorDouble getAnglePenalty(VectorDouble u_plus, VectorDouble u_minus);

}