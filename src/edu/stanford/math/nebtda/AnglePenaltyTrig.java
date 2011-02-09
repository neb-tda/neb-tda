package edu.stanford.math.nebtda;

/**
 * Computes an angle penalty on a 1-cell.
 * 
 * @author Henry Adams
 * @author Atanas Atanasov
 */
public class AnglePenaltyTrig extends AnglePenalty {

	protected double fnZero;
	protected double fnOne;
	
	/**
	 * AnglePenaltyTrig constructor.
	 * 
	 * @param fnZero				angle where angle forces first become nonzero
	 * @param fnOne					angle where angle forces first reach their maximum
	 */
	public AnglePenaltyTrig(double fnZero, double fnOne) {
		setFnZeroFnOne(fnZero, fnOne);
	}
	
	/**
	 * Sets the angle where angle forces first become nonzero, and the angle where angle forces first reach their maximum.
	 * 
	 * @param fnZero				angle where angle forces first become nonzero
	 * @param fnOne					angle where angle forces first reach their maximum
	 */
	public void setFnZeroFnOne(double fnZero, double fnOne) {
		if ((fnZero < 0) || (fnOne < 0)) {
			throw new IllegalArgumentException("The fnZero and fnOne constants must be nonnegative.");
		}
		if ((fnZero > Math.PI) || (fnOne > Math.PI)) {
			throw new IllegalArgumentException("The fnZero and fnOne constants must be less than or equal to pi.");
		}
		if (fnZero >= fnOne ) {
			throw new IllegalArgumentException("The constants must satisfy that fnZero is less than fnOne.");
		}
		this.fnZero = fnZero;
		this.fnOne = fnOne;
	}

	/**
	 * Gets the angle where angle forces first become nonzero.
	 * 
	 * @return fnZero
	 */
	public double getFnZero() {
		return fnZero;
	}
	
	/**
	 * Gets the angle where angle forces first reach their maximum.
	 * 
	 * @return fnOne
	 */
	public double getFnOne() {
		return fnOne;
	}

	/**
	 * Computes the angle penalty at a node.
	 * The angle penalty is h(theta)(u_plus - u_minus), where theta is the angle between  u_plus and u_minus.
	 * Function h is defined by h(x) = 0 if x <= fnZero
	 * 							h(x) = (1-cos(pi(x-fnZero)/(x-fnOne)))/2 if fnZero < x < fnOne
	 * 							h(x) = 1 if x >= fnOne. 
	 * 
	 * @param u_plus		first adjacent edge
	 * @param u_minus		second adjacent edge
	 * @return 				angle penalty
	 */
	public VectorDouble getAnglePenalty(VectorDouble u_plus, VectorDouble u_minus) {
		double angle = VectorDouble.angle(u_plus, u_minus);
		if (angle <= fnZero) {
			return new VectorDouble(u_plus.coordinates.length);
		} else if (angle < fnOne) {
			double penalty = (1 - Math.cos(Math.PI / (fnOne - fnZero) * (angle - fnZero))) / 2;
			// TODO: simplify expression, as in paper.
			return VectorDouble.linearCombination(penalty, u_plus, -penalty, u_minus);
		} else {
			return VectorDouble.subtract(u_plus, u_minus);
		}	
	}

}