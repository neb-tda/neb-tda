package edu.stanford.math.nebtda;

/**
 * Abstract class for generating initial bands.
 * 
 * @author Henry Adams
 * @author Atanas Atanasov
 */
public abstract class InitialBand {

	protected int numBandPoints;
	
	/**
	 * Constructor.
	 * 
	 * @param numBandPoints		number of points in the band
	 */
	public InitialBand(int numBandPoints) {
		setNumBandPoints(numBandPoints);
	}
	
	/**
	 * Sets the number of points in the band.
	 * 
	 * @param numBandPoints		number of points in the band
	 */
	public void setNumBandPoints(int numBandPoints) {		
		if (numBandPoints < 3) {
			throw new IllegalArgumentException("There must be at least 3 points in the band.");
		}
		this.numBandPoints = numBandPoints;
	}
	
	/**
	 * Gets the number of points in the band.
	 * 
	 * @return number of points in the band
	 */
	public int getNumBandPoints() {
		return numBandPoints;
	}
	
	/**
	 * Builds an initial band.
	 * 
	 * @param start		first endpoint
	 * @param end		second endpoint
	 * @return 			initial band
	 */
	public abstract VectorDouble[] getInitialBand(VectorDouble start, VectorDouble end);
	
	/**
	 * Builds a line segment.
	 * 
	 * @param start		first endpoint
	 * @param end		second endpoint
	 * @return 			line segment
	 */
	public VectorDouble[] line(VectorDouble start, VectorDouble end) {
		VectorDouble[] pts = new VectorDouble[numBandPoints];
		for (int i = 0; i < numBandPoints; i++) {
			pts[i] = VectorDouble.linearCombination(1.0 - ((double) i) / ((double) (numBandPoints - 1)), start,
					((double) i) / ((double) (numBandPoints - 1)), end);
		}
		return pts;
	}

	/**
	 * Builds a wedge, which is the union of two line segments defined by two endpoints and a corner.
	 * 
	 * @param start		first endpoint
	 * @param end		second endpoint
	 * @param corner	corner
	 * @return 			wedge
	 */
	public VectorDouble[] wedge(VectorDouble start, VectorDouble end, VectorDouble corner) {
		VectorDouble a = VectorDouble.subtract(corner, start);
		VectorDouble b = VectorDouble.subtract(corner, end);
		double normA = a.norm();
		double normB = b.norm();
		if ((normA == 0) || (normB == 0)) {
			throw new IllegalArgumentException("Input vectors are equal.");
		}
		double unit = (normA + normB) / (numBandPoints - 1);
		int indCutoff = (int) Math.floor(normA / unit);

		VectorDouble[] points = new VectorDouble[numBandPoints];
		points[0] = start;
		for (int i = 1; i <= indCutoff; i++) {
			points[i] = VectorDouble.linearCombination(1.0, points[i - 1], unit / normA, a);
		}
		points[numBandPoints - 1] = end;
		for (int i = numBandPoints - 2; i > indCutoff; i--) {
			points[i] = VectorDouble.linearCombination(1.0, points[i + 1], unit / normB, b);
		}
		return points;
	}

	/**
	 * Builds a circular arc which lies on the unique circle defined by VectorDoubles start, end, and dir.
	 * The arc passes through VectorDouble dir according to boolean passThroughDir.
	 * 
	 * @param start				first endpoint
	 * @param end				second endpoint
	 * @param dir				direction vector
	 * @param passThroughDir	if true, then the band passes through dir. If false, then the band passes through -dir.
	 * @return					circular arc
	 */
	public VectorDouble[] circleArc(VectorDouble start, VectorDouble end, VectorDouble dir, Boolean passThroughDir) {
		VectorDouble center = circumCenter(start, end, dir);
		VectorDouble shiftedStart = VectorDouble.subtract(start, center);
		VectorDouble shiftedEnd = VectorDouble.subtract(end, center);
		VectorDouble shiftedDir = VectorDouble.subtract(dir, center);
		VectorDouble perp = VectorDouble.perpendicularProjection(shiftedStart, shiftedEnd);
		if (perp.norm() == 0) {
			 throw new IllegalArgumentException("The norm of perp is zero.");
		}
		perp.scalarMult(shiftedStart.norm() / perp.norm());
		double angle = VectorDouble.angle(shiftedStart, shiftedEnd);
		 
		// Let alpha be the angle between shiftedStart - shiftedDir and shiftedEnd - shiftedDir. The sign of cosAlpha 
		// depends only on the dot product of shiftedStart - shiftedDir and shiftedEnd - shiftedDir. In the case
		// cosAlpha < 0, then alpha > PI/2, and shiftedDir is on the short path between shiftedStart and shiftedEnd.
		Boolean dirOnShortPath = (VectorDouble.dotProduct(VectorDouble.subtract(shiftedStart, shiftedDir), VectorDouble
				.subtract(shiftedEnd, shiftedDir)) < 0);
		if ((dirOnShortPath & (!passThroughDir)) || ((!dirOnShortPath) & passThroughDir)) {
			angle = 2 * Math.PI - angle;
			perp.negate();
		}
		
		VectorDouble[] points = new VectorDouble[numBandPoints];
		for (int i = 0; i < numBandPoints; i++) {
			points[i] = VectorDouble.linearCombination(Math.cos(angle * i / ((double) (numBandPoints - 1))), shiftedStart,
					Math.sin(angle * i / ((double) (numBandPoints - 1))), perp, 1.0, center);
		}
		return points;
		
	}
	
	/**
	 * First, builds a circular arc which lies on the unique circle defined by the unit vectors in the directions
	 * of VectorDoubles start, end, and dir. This arc passes through VectorDouble dir according to boolean 
	 * passThroughDir. Then, the points in this arc are scaled so that the norms of the points change linearly 
	 * from the norm of VectorDouble start to the norm of VectorDouble end.
	 * 
	 * @param start				first endpoint
	 * @param end				second endpoint
	 * @param dir				direction vector
	 * @param passThroughDir	if true, then the band passes through dir. If false, then the band passes through -dir.
	 * @return					circular arc
	 */
	public VectorDouble[] approxCircleArc(VectorDouble start, VectorDouble end, VectorDouble dir, Boolean passThroughDir) {
		VectorDouble normalizedStart = VectorDouble.unitVector(start);
		VectorDouble normalizedEnd = VectorDouble.unitVector(end);
		VectorDouble normalizedDir = VectorDouble.unitVector(dir);
		VectorDouble[] circleArc = circleArc(normalizedStart, normalizedEnd, normalizedDir, passThroughDir);
		
		for (int i = 0; i < circleArc.length; i++) {
			circleArc[i].scalarMult((1.0 - (double) i / (double) (circleArc.length - 1)) * start.norm() + ((double) i / (double) (circleArc.length - 1)) * end.norm());
		}
		return circleArc;
	}

	/**
	 * Returns the circumcenter of the three vectors.
	 * 
	 * @param a			first vector
	 * @param b			second vector
	 * @param c			third vector
	 * @return			the circumcenter
	 */
	public static VectorDouble circumCenter(VectorDouble a, VectorDouble b, VectorDouble c) {
		if ((a.coordinates.length != a.coordinates.length) || (b.coordinates.length != c.coordinates.length)) {
			throw new IllegalArgumentException("Input vectors do not have the same length.");
		}
		double r = VectorDouble.subtract(a, b).norm();
		double s = VectorDouble.subtract(b, c).norm();
		double t = VectorDouble.subtract(c, a).norm();
		if ((s == 0) || (t == 0) || (r == 0)) {
			throw new IllegalArgumentException("Input vectors are equal");
		}

		double barCoord1 = s * s * (-s * s + t * t + r * r);
		double barCoord2 = t * t * (s * s - t * t + r * r);
		double barCoord3 = r * r * (s * s + t * t - r * r);
		if (barCoord1 + barCoord2 + barCoord3 == 0) {
			throw new IllegalArgumentException("Input vectors are colinear");
		}
		VectorDouble center = VectorDouble.linearCombination(barCoord1, a, barCoord2, b, barCoord3, c);
		center.scalarMult(1.0 / (barCoord1 + barCoord2 + barCoord3));
		return center;
	}

}
