package edu.stanford.math.nebtda;

/**
 * Implements a node in a higher dimensional cell.
 * 
 * @author Henry Adams
 * @author Atanas Atanasov
 */
public class CellHigherNode {
	
	protected static int count = 1;
	
	protected VectorDouble point;
	protected VectorDouble force;
	protected boolean isBoundary;
	protected int id; // identification number of the node.

	/**
	 * Constructor.
	 * 
	 * @param point			location of the node
	 * @param isBoundary	true if the node is a boundary node
	 */
	public CellHigherNode(VectorDouble point, boolean isBoundary) {
		if (point == null) {
			throw new IllegalArgumentException("Point reference is null.");
		}
		setPoint(point);
		this.force = null;
		this.isBoundary = isBoundary;
		id = (count++);
	}
	
	/**
	 * Sets the point.
	 * 
	 * @param point		location of the node
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
	 * @return location of the node
	 */
	public VectorDouble getPoint() {
		return point;
	}
	
	/**
	 * Returns the id number of the node as a string.
	 * 
	 * @return id
	 */
	public String toString() {
		return Integer.toString(id);
	}
	
	/**
	 * Sets the count index.
	 * 
	 * @param count		the new count index
	 */
	public void setCount(int count) {
		CellHigherNode.count = count;
	}
	
}
