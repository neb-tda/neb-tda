package edu.stanford.math.nebtda;

// TODO: made some changes to movement around loop. Not yet fully tested. Test!
/**
 * Implements the boundary of a 2-cell as a piecewise linear loop.
 * 
 * @author Henry Adams
 * @author Atanas Atanasov
 */
public class BoundaryLoop extends Boundary {

	protected VectorDouble[] points;

	/**
	 * Constructor.
	 * 
	 * @param points				points defining the loop
	 */
	public BoundaryLoop(VectorDouble[] points) {
		setPoints(points);
	}

	/**
	 * Constructor.
	 * 
	 * @param bands					bands around the loop
	 * @param orientations			orientations of the bands; true is forward
	 */
	public BoundaryLoop(VectorDouble[][] bands, boolean[] orientations) {
		this(buildPoints(bands, orientations));
	}

	/**
	 * Constructor.
	 * 
	 * @param bands					bands around the loop
	 * @param orientations			orientations of the bands; true is forward
	 */
	public BoundaryLoop(Cell1[] bands, boolean[] orientations) {
		this(buildPoints(bands, orientations));
	}

	/**
	 * Sets the points defining the loop.
	 * 
	 * @param points	points defining the loop
	 */
	public void setPoints(VectorDouble[] points) {
		VectorDouble.checkConsistentDimension(points);
		if (points.length < 3) {
			throw new IllegalArgumentException("There are less than 3 points.");
		}
		
		for (int i = 0; i < points.length; i++) {
			if (VectorDouble.equals(points[i], points[(i + 1) % (points.length)])) {
				throw new IllegalArgumentException("Two consecutive points on the loop are equal.");
			}
		}	
		this.points = points;
	}

	/**
	 * Gets the points defining the loop.
	 * 
	 * @return points defining the loop
	 */
	public VectorDouble[] getPoints() {
		return points;
	}

	/**
	 * Builds the loop points.
	 * 
	 * @param bands				bands defining the loop
	 * @param orientations		orientations of the bands
	 * @return the loop points
	 */
	private static VectorDouble[] buildPoints(VectorDouble[][] bands, boolean[] orientations) {
		int i, j, length;
		VectorDouble end, start;

		if (bands == null) {
			throw new IllegalArgumentException("The bands reference is null.");
		}
		if (orientations == null) {
			throw new IllegalArgumentException("The orientations reference is null.");
		}
		if (bands.length != orientations.length) {
			throw new IllegalArgumentException("The bands and orientations are inconsistent.");
		}

		length = 0;
		for (i = 0; i < bands.length; i++) {
			if (bands[i] == null) {
				throw new IllegalArgumentException("Some band reference is null.");
			}
			if (bands[i].length < 2) {
				throw new IllegalArgumentException("All bands should have at least 2 nodes.");
			}
			length += bands[i].length;
		}
		length -= bands.length;

		for (i = 0; i < bands.length; i++) {
			j = (i + 1) % (bands.length);

			end = (orientations[i]) ? (bands[i][bands[i].length - 1]) : (bands[i][0]);
			start = (orientations[j]) ? (bands[j][0]) : (bands[j][bands[j].length - 1]);
			
			if (!VectorDouble.equals(end, start)) {
				throw new IllegalArgumentException("The given bands do not glue together.");
			}
		}

		VectorDouble[] points = new VectorDouble[length];
		int count = 0;

		for (i = 0; i < bands.length; i++) {
			if (orientations[i]) {
				for (j = 0; j < bands[i].length - 1; j++, count++) {
					points[count] = new VectorDouble(bands[i][j]);
				}
			} else {
				for (j = bands[i].length - 1; j > 0; j--, count++) {
					points[count] = new VectorDouble(bands[i][j]);
				}
			}
		}
		
		return points;
	}

	/**
	 * Builds the loop points.
	 * 
	 * @param bands				bands defining the loop
	 * @param orientations		orientations of the bands
	 * @return the loop points
	 */
	public static VectorDouble[] buildPoints(Cell1[] bands, boolean[] orientations) {
		int i, j;
		VectorDouble[][] bands_points;

		if (bands == null) {
			throw new IllegalArgumentException("The bands reference is null.");
		}

		bands_points = new VectorDouble[bands.length][];

		for (i = 0; i < bands.length; i++) {
			if (bands[i] == null) {
				throw new IllegalArgumentException("Some band reference is null.");
			}

			bands_points[i] = new VectorDouble[bands[i].points.length];

			for (j = 0; j < bands_points[i].length; j++) {
				bands_points[i][j] = bands[i].points[j];
			}
		}

		return buildPoints(bands_points, orientations);
	}

	/**
	 * Moves a point along a geodesic in the boundary. If the point is not in the 
	 * boundary loop, then it is first moved to the nearest location in the loop.
	 * 
	 * @param point			initial point on the boundary
	 * @param tangent		initial direction to move
	 * @return 				result after moving along a geodesic in the boundary
	 */
	public VectorDouble followGeodesic(VectorDouble point, VectorDouble tangent) {
		VectorDouble tmp, dir;
		double usage = 0.0;
		int segment;
		boolean forward;

		segment = movePointAndGetSegment(point);

		for (;;) {
			dir = VectorDouble.subtract(points[(segment + 1) % (points.length)], points[segment]);
			tmp = VectorDouble.projection(dir, tangent);
			usage = VectorDouble.dotProduct(tmp, dir);

			if (usage == 0.0) {
				return point;
			} else {
				if (usage > 0.0) {
					dir = VectorDouble.subtract(points[(segment + 1) % (points.length)], point);
					forward = true;
				} else {
					dir = VectorDouble.subtract(points[segment], point);
					forward = false;
				}
			}

			usage = dir.normSquared();

			if (usage == 0.0) {
				return point;
			}

			usage = VectorDouble.dotProduct(dir, tmp) / usage;

			if (usage <= 1.0) {
				point.add(tmp);
				return point;
			} else {
				tmp.scalarMult((usage - 1) / usage);

				if (forward) {
					segment = (segment + 1) % (points.length);
					point.assign(points[segment]);
				} else {
					point.assign(points[segment]);
					segment = (segment + points.length - 1) % (points.length);
				}
			}
		}
	}

	/**
	 * Projects a vector to the tangent space of the boundary. If the base point is not in the 
	 * boundary loop, then it is first moved to the nearest location in the loop.
	 * 
	 * @param point			point on the boundary
	 * @param vector		vector
	 * @return 				vector projected to the tangent space of the boundary at the point
	 */
	public VectorDouble projectTangent(VectorDouble point, VectorDouble vector) {
		int segment = movePointAndGetSegment(point);
		return VectorDouble.projection(VectorDouble.subtract(points[segment], points[(segment + 1) % (points.length)]),
				vector);
	}
	
	/**
	 * Moves the point to the nearest location in the loop, and gets the segment index of that location.
	 * 
	 * @param point			point on the boundary
	 * @return 				segment index
	 */
	public int movePointAndGetSegment(VectorDouble point) {
		int segment = 0;
		double min = VectorDouble.distanceToSegment(point, points[0], points[1]);
		double distance;
		for (int i = 1; i < points.length; i++) {
			distance = VectorDouble.distanceToSegment(point, points[i], points[(i + 1) % (points.length)]);
			if (distance < min) {
				min = distance;
				segment = i;
			}
		}

		point.moveToSegment(points[segment], points[(segment + 1) % (points.length)]);
		return segment;
	}
	
	/**
	 * Computes the average of the points defining the loop.
	 * 
	 * @return the average of the points defining the loop
	 */
	public VectorDouble computeAvgCenter() {
		VectorDouble center = new VectorDouble(points[0]);
		for (int i = 1; i < points.length; i++) {
			center.add(points[i]);
		}
		center.scalarMult(1.0 / (double) points.length);
		return center;
	}

}
