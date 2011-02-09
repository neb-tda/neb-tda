package edu.stanford.math.nebtda;

import java.util.Random;

/**
 * Generates initial bands for general data sets in Euclidean space.
 * 
 * @author Henry Adams
 * @author Atanas Atanasov
 */
public class InitialBandEuclidean extends InitialBand {
	
	protected static Random r = new Random();

	/**
	 * Constructor.
	 * 
	 * @param numBandPoints		number of points in the band
	 */
	public InitialBandEuclidean(int numBandPoints) {
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
		int dim = start.coordinates.length;	
		VectorDouble midpoint = VectorDouble.linearCombination(0.5, start, 0.5, end);
		VectorDouble shiftedStart = VectorDouble.subtract(start, midpoint);
		VectorDouble dir = VectorDouble.perpendicularProjection(shiftedStart, DistributionNormal.nextSphereVector(dim));
		while (dir.norm() == 0) {
			dir = VectorDouble.perpendicularProjection(shiftedStart, DistributionNormal.nextSphereVector(dim));
		}
		
		dir.scalarMult(r.nextDouble() * shiftedStart.norm() / dir.norm());
		dir.add(midpoint);
		return circleArc(start, end, dir, true);

	}
	
}