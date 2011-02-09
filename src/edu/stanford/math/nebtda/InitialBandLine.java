package edu.stanford.math.nebtda;

/**
 * Generates a straight line initial band.
 * 
 * @author Henry Adams
 * @author Atanas Atanasov
 */
public class InitialBandLine extends InitialBand {

	/**
	 * Constructor.
	 * 
	 * @param numBandPoints		number of points in the band
	 */
	public InitialBandLine(int numBandPoints) {
		super(numBandPoints);
	}

	/**
	 * Builds a straight line initial band.
	 * 
	 * @param start		first endpoint
	 * @param end		second endpoint
	 * @return 			initial band
	 */
	@Override
	public VectorDouble[] getInitialBand(VectorDouble start, VectorDouble end) {
		return line(start, end);
	}
	
}
