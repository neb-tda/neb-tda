package edu.stanford.math.nebtda;

import java.io.PrintWriter;

/**
 * Implements a 0-cell.
 * 
 * @author Henry Adams
 * @author Atanas Atanasov
 */
public class Cell0 extends Cell implements MetricSpace<Cell0>{
		
	public static Cell0 instance;
	static {
		instance = new Cell0();
	}

	protected VectorDouble point;
	
	/**
	 * Empty constructor.
	 */
	private Cell0() {
	}
	
	/**
	 * Constructor.
	 * 
	 * @param point		location of the 0-cell 
	 */
	public Cell0(VectorDouble point) {
		setPoint(point);
	}
	
	/**
	 * Sets the point.
	 * 
	 * @param point		location of the 0-cell
	 */
	public void setPoint(VectorDouble point) {
		if (point == null) {
			throw new IllegalArgumentException("Point reference is null.");
		}
		this.point = point;
	}
	
	/**
	 * Gets the point.
	 * 
	 * @return location of the 0-cell
	 */
	public VectorDouble getPoint() {
		return point;
	}
	
	/**
	 * Computes the density value.
	 * 
	 * @param density	density estimator
	 * @return 			density value
	 */
	public double computeValue(Density density) {
		return density.getValue(point);
	}
	
	/**
	 * Runs the 0-cell using the mean shift method.
	 * 
	 * @param density		density estimator
	 * @param p0			trial parameters
	 * @param verbose		if true, prints trial to screen
	 * @param pointWriter	writes the point locations
	 * @param distanceWriter	writes the forces
	 */
	public void run(Density density, Parameters0 p0, boolean verbose, PrintWriter pointWriter, PrintWriter distanceWriter) {
		
		VectorDouble meanShift = density.getMeanShift(point);
		double distanceToMove = VectorDouble.distance(point, meanShift);
		int count = 0;
		
		if (verbose) {
			System.out.println("Step " + count);
			System.out.println(point.toString(false, false));
			System.out.println("Distance to move  = " + distanceToMove);
		}
		pointWriter.println(point.toString(false, false));
		distanceWriter.println(distanceToMove);
		
		while (distanceToMove > p0.convergence) {	
			count ++;
			point = meanShift;
			meanShift = density.getMeanShift(point);
			distanceToMove = VectorDouble.distance(point, meanShift);
			
			if (verbose) {
				System.out.println("Step " + count);
				System.out.println(point.toString(false, false));
				System.out.println("Distance to move  = " + distanceToMove);
			}
			pointWriter.println(point.toString(false, false));
			distanceWriter.println(distanceToMove);
		}
		setValue(density.getValue(point));
	}
	
	/**
	 * Computes the distance between two 0-cells.
	 * 
	 * @param a		first 0-cell
	 * @param b		second 0-cell
	 * @return 		distance
	 */
	public double distance(Cell0 a, Cell0 b) {
		return VectorDouble.subtract(a.point, b.point).norm();
	}
	
}
