package edu.stanford.math.nebtda;

import java.io.PrintWriter;

/**
 * Implements a 1-cell.
 * 
 * @author Henry Adams
 * @author Atanas Atanasov
 */
public class Cell1 extends Cell implements MetricSpace<Cell1>{
	
	public static Cell1 instance;
	static {
		instance = new Cell1();
	}
	protected static Tangent tangentInstance = new TangentNaive();

	protected VectorDouble[] points;
	protected VectorDouble[] forces;
	protected boolean converged;
	
	/**
	 * Empty constructor.
	 */
	public Cell1() {
	}
	
	/**
	 * Constructor.
	 * 
	 * @param points	location of the 1-cell 
	 */
	public Cell1(VectorDouble[] points) {
		setPoints(points);
		forces = new VectorDouble[points.length];
		converged = false;
	}
	
	/**
	 * Sets the points.
	 * 
	 * @param points	location of the 1-cell
	 */
	public void setPoints(VectorDouble[] points) {
		VectorDouble.checkConsistentDimension(points);
		this.points = points;
	}
	
	/**
	 * Gets the points.
	 * 
	 * @return location of the 1-cell
	 */
	public VectorDouble[] getPoints() {
		return points;
	}
	
	/**
	 * Computes the density value.
	 * 
	 * @param density	density estimator
	 * @return 			density value
	 */
	public double computeValue(Density density) {
		double result = density.getValue(points[0]);
		double cur;
		for (int i = 1; i < points.length; i++) {
			cur = density.getValue(points[i]);
			if (cur < result) {
				result = cur;
			}
		}
		return result;
	}
	
	/**
	 * Runs the 1-cell. If the cell converges, the density value is then calculated and boolean converged is set to true.
	 * 
	 * @param density		density estimator
	 * @param p1			trial parameters
	 * @param verbose		if true, prints trial to screen
	 * @param pointWriter	writes the point locations
	 * @param forceWriter	writes the forces
	 */
	public void run(Density density, Parameters1 p1, boolean verbose, PrintWriter pointWriter, PrintWriter forceWriter) {
		int i, count = 0;
		
		if (verbose) { System.out.println("Step " + count); }
		for (i = 0; i < points.length; i++) {
			if (verbose) { System.out.println(points[i].toString(false, false)); }
			pointWriter.println(points[i].toString(false, false));
		}

		double avgForceNorms = computeForces(density, p1, verbose, forceWriter);
		while (avgForceNorms > p1.convergence && count < p1.maxSteps) {
			for (i = 1; i < points.length - 1; i++) {
				points[i].add(VectorDouble.scalarMult(p1.step, forces[i]));
			}
			
			count++;
			if (verbose) {System.out.println("Step " + count); }
			for (i = 0; i < points.length; i++) {
				if (verbose) { System.out.println(points[i].toString(false, false)); }
				pointWriter.println(points[i].toString(false, false));
			}
			
			avgForceNorms = computeForces(density, p1, verbose, forceWriter);
		}
		if (avgForceNorms <= p1.convergence) {
			setValue(computeValue(density));
			converged = true;
		}
	}
	
	/**
	 * Computes the forces at each node.
	 * 
	 * @param density		density estimator
	 * @param p1			trial parameters
	 * @param verbose		if true, prints trial to screen
	 * @param forceWriter	writes the forces
	 * @return 				the average force norm
	 */
	public double computeForces(Density density, Parameters1 p1, boolean verbose, PrintWriter forceWriter) {
		VectorDouble tangent, u_plus = null, u_minus = null;
		VectorDouble gradF;
		VectorDouble springF;
		VectorDouble angleF;
		double gradFNorms = 0.0;
		double springFNorms = 0.0;
		double angleFNorms = 0.0;
		double totalFNorms = 0.0;
		
		for (int i = 1; i < points.length - 1; i++) {
			tangent = tangentInstance.getTangent(points[i - 1], points[i], points[i + 1]);
			u_plus = VectorDouble.subtract(points[i + 1], points[i]);
			u_minus = VectorDouble.subtract(points[i], points[i - 1]);

			gradF = VectorDouble.perpendicularProjection(tangent, density.getGradient(points[i]));
			gradF.scalarMult(p1.gradConst);
			gradFNorms += gradF.norm();

			springF = VectorDouble.scalarMult((u_plus.norm() - u_minus.norm()) / tangent.norm(), tangent);
			springFNorms += springF.norm();

			angleF = p1.anglePenalty.getAnglePenalty(u_plus, u_minus);
			angleFNorms += angleF.norm();
			forces[i] = VectorDouble.linearCombination(1.0, gradF, 1.0, springF, 1.0, angleF);
			totalFNorms += forces[i].norm();
		}

		double averageGradFNorm = gradFNorms / (double) (points.length - 2);
		double averageSpringFNorm = springFNorms / (double) (points.length - 2);
		double averageAngleFNorm = angleFNorms / (double) (points.length - 2);
		double averageTotalFNorm = totalFNorms / (double) (points.length - 2);
		
		if (verbose) {
			System.out.println("Average gradient force = " + averageGradFNorm);
			System.out.println("Average spring force   = " + averageSpringFNorm);
			System.out.println("Average angle force    = " + averageAngleFNorm);
			System.out.println("Average total force    = " + averageTotalFNorm);
		}
		forceWriter.print(averageGradFNorm + ", ");
		forceWriter.print(averageSpringFNorm + ", ");
		forceWriter.print(averageAngleFNorm + ", ");
		forceWriter.println(averageTotalFNorm);

		return averageTotalFNorm;
	}

	/**
	 * Returns true if the band is within the nearness threshold of the vertices.
	 * 
	 * @param vertices				vertices
	 * @param nearnessThreshold		nearness threshold
	 * @return						true or false
	 */
	public boolean bandNearVertices(VectorDouble[] vertices, double nearnessThreshold) {

		for (int i = 0; i < vertices.length; i++) {
			if (distanceToBand(vertices[i]) < nearnessThreshold) {
				return true;
			}
		}
		return false;
	}
	
	/**
	 * Returns the distance from the vertex to the band.
	 * 
	 * @param vertex	vertex
	 * @return 		distance to the vertex
	 */
	public double distanceToBand(VectorDouble vertex) {
		double cur;
		double distance = VectorDouble.distanceToSegment(vertex, points[0], points[1]);
		for (int i = 1; i < points.length - 1; i++) {
			cur = VectorDouble.distanceToSegment(vertex, points[i], points[i + 1]);
			if (cur < distance) {
				distance = cur;
			}
		}
		return distance;
	}
	
	/**
	 * Computes the distance between two 1-cells.
	 * 
	 * @param a		first 1-cell
	 * @param b		second 1-cell
	 * @return 		distance
	 */
	public double distance(Cell1 a, Cell1 b) {
		int i, size;
		size = b.points.length;

		if (a.points.length != size) {
			throw new IllegalArgumentException("The band lengths are inconsistent.");
		}

		double distanceFirst, distanceLast, distancesSum;
		distanceFirst = VectorDouble.subtract(a.points[0], b.points[0]).norm();
		distanceLast = VectorDouble.subtract(a.points[size - 1], b.points[size - 1]).norm();
		distancesSum = 0;

		for (i = 1; i < size - 1; i++) {
			distancesSum += VectorDouble.subtract(a.points[i], b.points[i]).norm();
		}

		if ((distanceFirst > 0) || (distanceLast > 0)) {
			return (distanceFirst + distanceLast + distancesSum) / size;
		} else {
			return distancesSum / (size - 2);
		}
	}
	
}
