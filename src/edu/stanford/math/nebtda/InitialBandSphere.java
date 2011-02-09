package edu.stanford.math.nebtda;

/**
 * Generates initial bands for datasets which are normalized to lie on a unit sphere.
 * 
 * @author Henry Adams
 * @author Atanas Atanasov
 */
public class InitialBandSphere extends InitialBand {

	/**
	 * Constructor.
	 * 
	 * @param numBandPoints		number of points in the band
	 */
	public InitialBandSphere(int numBandPoints) {
		super(numBandPoints);
	}

	/**
	 * Builds an initial band.
	 * 
	 * @param start		first endpoint
	 * @param end		second endpoint
	 * @return 			initial band
	 */
	@Override
	public VectorDouble[] getInitialBand(VectorDouble start, VectorDouble end) {
		return getInitialBand(start, end, DistributionNormal.nextSphereVector(start.getDimension()));
	}
	
	/**
	 * Builds an initial band.
	 * 
	 * @param start		first endpoint
	 * @param end		second endpoint
	 * @param dir		specifies a direction
	 * @return 			initial band
	 */
	public VectorDouble[] getInitialBand(VectorDouble start, VectorDouble end, VectorDouble dir) {
		return approxCircleArc(start, end, dir, false);
	}
	
}
