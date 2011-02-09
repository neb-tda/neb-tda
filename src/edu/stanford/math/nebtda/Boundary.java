package edu.stanford.math.nebtda;

/**
 * Abstract class for the boundary of a higher dimensional cell.
 * 
 * @author Henry Adams
 * @author Atanas Atanasov
 */
public abstract class Boundary {

	/**
	 * Moves a point along a geodesic in the boundary.
	 * 
	 * @param point			initial point on the boundary
	 * @param tangent		initial direction to move
	 * @return 				result after moving along a geodesic in the boundary
	 */
	public abstract VectorDouble followGeodesic(VectorDouble point, VectorDouble tangent);

	/**
	 * Projects a vector to the tangent space of the boundary.
	 * 
	 * @param point			point on the boundary
	 * @param vector		vector
	 * @return 				vector projected to the tangent space of the boundary at the point
	 */
	public abstract VectorDouble projectTangent(VectorDouble point, VectorDouble vector);
	
}
