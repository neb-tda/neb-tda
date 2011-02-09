package edu.stanford.math.nebtda;

import java.io.PrintWriter;
import java.util.Iterator;
import java.util.Random;
import java.util.Set;

import org.jgrapht.DirectedGraph;
import org.jgrapht.UndirectedGraph;
import org.jgrapht.graph.AsUndirectedGraph;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultEdge;

/**
 * Implements a higher dimensional cell.
 * 
 * @author Henry Adams
 * @author Atanas Atanasov
 */
public class CellHigher extends Cell {

	protected static Random r = new Random();

	protected UndirectedGraph<CellHigherNode, DefaultEdge> graph;
	protected int cellDimension;
	protected Boundary boundary;
	protected int numMovingNodes;
	boolean converged;

	/**
	 * Constructor. If boundary is null, then the boundary nodes are fixed.
	 * 
	 * @param graph				graph representing the cell
	 * @param cellDimension		cell dimension
	 * @param boundary 			boundary
	 */
	public CellHigher(UndirectedGraph<CellHigherNode, DefaultEdge> graph, int cellDimension, Boundary boundary) {
		setGraph(graph);
		setMeshDimension(cellDimension);
		this.boundary = boundary;
		this.converged = false;
	}
	
	/**
	 * Sets the graph and the number of moving nodes.
	 * 
	 * @param graph				graph representing the cell
	 */
	public void setGraph(UndirectedGraph<CellHigherNode, DefaultEdge> graph) {
		if (graph == null) {
			throw new IllegalArgumentException("The graph reference is null.");
		}
		this.graph = graph;
		numMovingNodes = 0;
		for (CellHigherNode node : graph.vertexSet()) {
			if (!node.isBoundary || boundary != null) {
				numMovingNodes++;
			}
		}
	}
	
	/**
	 * Gets the graph.
	 * 
	 * @return graph representing the cell
	 */
	public UndirectedGraph<CellHigherNode, DefaultEdge> getGraph() {
		return graph;
	}

	/**
	 * Sets the mesh dimension.
	 * 
	 * @param cellDimension		cell dimension
	 */
	public void setMeshDimension(int cellDimension) {
		if (cellDimension <= 0) {
			throw new IllegalArgumentException("Mesh dimension is not positive.");
		}
		this.cellDimension = cellDimension;
	}

	/**
	 * Gets the cell dimension.
	 * 
	 * @return cell dimension
	 */
	public int getMeshDimension() {
		return cellDimension;
	}
	
	/**
	 * Computes the density value.
	 * 
	 * @param density	density estimator
	 * @return 			density value
	 */
	public double computeValue(Density density) {
		Iterator<CellHigherNode> it = graph.vertexSet().iterator();
		if (!it.hasNext()) {
			throw new IllegalArgumentException("The cell is empty.");
		}
		double result = density.getValue(it.next().point);
		double cur;
		while (it.hasNext()) {
			cur = density.getValue(it.next().point);
			if (cur < result) {
				result = cur;
			}
		}
		return result;
	}
	
	/**
	 * Runs the cell. If the cell converges, the density value is then calculated and boolean converged is set to true.
	 * 
	 * @param density		density estimator
	 * @param p				trial parameters
	 * @param verbose		if true, prints trial to screen
	 * @param pointWriter	writes the point locations
	 * @param forceWriter	writes the forces
	 */
	public void run(Density density, ParametersHigher p, boolean verbose, PrintWriter pointWriter, PrintWriter forceWriter) {
		int count = 0;
		
		if (verbose) {
			System.out.println("Step " + count);
		}
		for (CellHigherNode node : graph.vertexSet()) {
			pointWriter.println(node.point.toString(false, false));
		}
		
		double avgForceNorms = computeForces(density, p, verbose, forceWriter);
		while (avgForceNorms > p.convergence && count < p.maxSteps) {
			for (CellHigherNode node : graph.vertexSet()) {
				if ((!node.isBoundary) || (boundary != null)) {
					if (node.isBoundary) { 
						node.point = boundary.followGeodesic(node.point, VectorDouble.scalarMult(p.step, node.force));
					} else {
						node.point.add(VectorDouble.scalarMult(p.step, node.force));
					}
					pointWriter.println(node.point.toString(false, false));
				}
			}
			
			count ++;
			if (verbose) {
				System.out.println("Step " + count);
			}
			
			avgForceNorms = computeForces(density, p, verbose, forceWriter);
		}
		if (avgForceNorms <= p.convergence) {
			setValue(computeValue(density));
			converged = true;
		}
	}

	/**
	 * Computes the forces at each node.
	 * 
	 * @param density		density estimator
	 * @param p				trial parameters
	 * @param verbose		if true, prints trial to screen
	 * @param forceWriter	writes the forces
	 * @return 				the average force norm
	 */
	public double computeForces(Density density, ParametersHigher p, boolean verbose, PrintWriter forceWriter) {
		VectorDouble point, gradF, springF;
		VectorDouble[] incidentPoints, tangentSpace;
		double sumGradFNorms = 0.0;
		double sumSpringFNorms = 0.0;
		double sumTotalFNorms = 0.0;

		for (CellHigherNode node : graph.vertexSet()) {
			if ((!node.isBoundary) || (boundary != null)) {
				incidentPoints = getIncidentNodes(node);
				point = node.point;
				if (node.isBoundary) { 
					gradF = VectorDouble.scalarMult(p.gradConst, boundary.projectTangent(point, density.getGradient(point)));
				} else {
					tangentSpace = VectorDouble.getTangentSpace(point, incidentPoints, cellDimension);
					gradF = VectorDouble.scalarMult(p.gradConst, VectorDouble.perpendicularProjection(tangentSpace, density.getGradient(point), true));
				}
				springF = getSpringForce(point, incidentPoints);
				sumGradFNorms += gradF.norm();
				sumSpringFNorms += springF.norm();
				node.force = VectorDouble.add(gradF, springF);
				sumTotalFNorms += node.force.norm();
			}
		}
		
		double averageGradFNorm = sumGradFNorms / numMovingNodes;
		double averageSpringFNorm = sumSpringFNorms / numMovingNodes;
		double averageTotalFNorm = sumTotalFNorms / numMovingNodes;
		
		if (verbose) {
			System.out.println("Average gradient force = " + averageGradFNorm);
			System.out.println("Average spring force   = " + averageSpringFNorm);
			System.out.println("Average total force    = " + averageTotalFNorm);
		}
		forceWriter.print(averageGradFNorm + ", ");
		forceWriter.print(averageSpringFNorm + ", ");
		forceWriter.println(averageTotalFNorm);

		return averageTotalFNorm;
	}
	
	/**
	 * Finds the incident nodes in the graph.
	 * 
	 * @param node		the given node
	 * @return			the incident nodes
	 */
	public VectorDouble[] getIncidentNodes(CellHigherNode node) {
		CellHigherNode source, target;
		
		Set<DefaultEdge> edges = graph.edgesOf(node);
		if (edges.size() == 0) {
			throw new IllegalArgumentException("The graph should have no isolated nodes.");
		}
		VectorDouble[] incidentNodes = new VectorDouble[edges.size()];
		
		int count = 0;
		for (DefaultEdge edge : edges) {
			source = graph.getEdgeSource(edge);
			target = graph.getEdgeTarget(edge);
			if (source == node) {
				incidentNodes[count] = target.point;
			} else {
				if (target == node) {
					incidentNodes[count] = source.point;
				} else {
					throw new IllegalArgumentException("The structure graph is inconsistent.");
				}
			}
			count++;
		}

		return incidentNodes;
	}
	
	/**
	 * Gets the spring force at a node.
	 * 
	 * @param point				the position of the node
	 * @param incidentPoints	the positions of the incident nodes
	 * @return					the spring force
	 */
	public VectorDouble getSpringForce(VectorDouble point, VectorDouble[] incidentPoints) {
		VectorDouble result = new VectorDouble(point.coordinates.length);
		for (int i = 0; i < incidentPoints.length; i++) {
			result.add(VectorDouble.subtract(incidentPoints[i], point));
		}
		return result;
	}

	/**
	 * Computes evenly spaced points about a loop.
	 * 
	 * @param vertices		vertices defining the loop
	 * @param numPoints		number of points
	 * @param randomize		if true, the position of the initial point is random
	 * @return				evenly spaced points about the loop
	 */
	private static VectorDouble[] buildAngularPoints(VectorDouble[] vertices, int numPoints, boolean randomize) {
		VectorDouble[] points;
		double step, tmp, advance;
		double[] lengths;
		int i, j;

		if (numPoints < 3) {
			throw new IllegalArgumentException("The number of points should be at least 3.");
		}

		lengths = new double[vertices.length];
		points = new VectorDouble[numPoints];
		step = 0.0;

		for (i = 0; i < vertices.length; i++) {
			lengths[i] = VectorDouble.distance(vertices[i], vertices[(i + 1) % (vertices.length)]);
			step += lengths[i];
		}

		step /= numPoints;
		tmp = (randomize) ? (r.nextDouble() * step) : (0.0);
		advance = 0.0;

		for (i = 0, j = 0; i < numPoints; i++) {
			while (lengths[j] <= tmp + advance) {
				tmp -= lengths[j] - advance;
				advance = 0.0;
				j = (j + 1) % (lengths.length);
			}

			advance += tmp;
			tmp = advance / lengths[j];
			points[i] = VectorDouble.linearCombination(1 - tmp, vertices[j], tmp, vertices[(j + 1) % (vertices.length)]);
			tmp = step;
		}

		return points;
	}

	/**
	 * Builds a web-shaped graph.
	 * 
	 * @param loop				boundary loop
	 * @param angularSize		number of points around the boundary
	 * @param radialSize		number of rings in the web
	 * @param center			location of the center of the web
	 * @param placeCenter		if true, the center is included in the web
	 * @return					web-shaped graph
	 */
	public static UndirectedGraph<CellHigherNode, DefaultEdge> buildWeb(BoundaryLoop loop, int angularSize, int radialSize,
			VectorDouble center, boolean placeCenter) {
		DirectedGraph<CellHigherNode, DefaultEdge> directed;
		CellHigherNode[][] nodes;
		CellHigherNode centerNode;
		VectorDouble[] points;
		VectorDouble point;
		double tmp;
		int i, j;

		if (loop == null) {
			throw new IllegalArgumentException("The loop reference is null.");
		}

		if (angularSize < 3) {
			throw new IllegalArgumentException("The angular size should be at least 3.");
		}

		if (radialSize < 1) {
			throw new IllegalArgumentException("The radial size should be at least 1.");
		}

		points = buildAngularPoints(loop.points, angularSize, false);
		nodes = new CellHigherNode[angularSize][radialSize];
		directed = new DefaultDirectedGraph<CellHigherNode, DefaultEdge>(DefaultEdge.class);

		for (i = 0; i < angularSize; i++) {
			for (j = 0; j < radialSize; j++) {
				tmp = ((double) (j + 1)) / ((double) radialSize);
				point = VectorDouble.linearCombination(tmp, points[i], 1 - tmp, center);
				nodes[i][j] = new CellHigherNode(point, (j == radialSize - 1));
				directed.addVertex(nodes[i][j]);
			}
		}

		for (i = 0; i < angularSize; i++) {
			for (j = 0; j < radialSize; j++) {
				directed.addEdge(nodes[i][j], nodes[(i + 1) % (angularSize)][j]);
				if (j < radialSize - 1) {
					directed.addEdge(nodes[i][j], nodes[i][j + 1]);
				}
			}
		}

		if (placeCenter) {
			centerNode = new CellHigherNode(new VectorDouble(center), false);
			directed.addVertex(centerNode);

			for (i = 0; i < angularSize; i++) {
				directed.addEdge(centerNode, nodes[i][0]);
			}
		}

		return new AsUndirectedGraph<CellHigherNode, DefaultEdge>(directed);
	}

	/**
	 * Builds a web-shaped graph. Cubic spline interpolation is used so that the web lies near the tangent space of the center.
	 * 
	 * @param loop				boundary loop
	 * @param angularSize		number of points around the boundary
	 * @param radialSize		number of rings in the web
	 * @param center			location of the center of the web
	 * @param placeCenter		if true, the center is included in the web
	 * @param lambda			scales the slant of the cone
	 * @param boundaryScale		scales the spline velocity at the boundary
	 * @param centerScale		scales the spline velocity at the center
	 * @return					web-shaped graph
	 */
	public static UndirectedGraph<CellHigherNode, DefaultEdge> buildSplineWeb(BoundaryLoop loop, int angularSize,
			int radialSize, VectorDouble center, boolean placeCenter, double lambda, double boundaryScale,
			double centerScale) {
		DirectedGraph<CellHigherNode, DefaultEdge> directed;
		CellHigherNode[][] nodes;
		CellHigherNode centerNode;
		VectorDouble[] vertices, points, boundaryTangent, centerTangent, tangentSpace;
		VectorDouble point, tip, centerMass;
		double t;
		int i, j;

		if (loop == null) {
			throw new IllegalArgumentException("The loop reference is null.");
		}

		if (angularSize < 3) {
			throw new IllegalArgumentException("The angular size should be at least 3.");
		}

		if (radialSize < 1) {
			throw new IllegalArgumentException("The radial size should be at least 1.");
		}

		if ((boundaryScale <= 0) || (centerScale <= 0)) {
			throw new IllegalArgumentException("Both scales should be positive.");
		}

		centerMass = new VectorDouble(center.coordinates.length);
		vertices = new VectorDouble[loop.points.length];

		for (i = 0; i < vertices.length; i++) {
			vertices[i] = loop.points[i];
			centerMass.linearCombination(1, 1 / ((double) vertices.length), vertices[i]);
		}

		tangentSpace = VectorDouble.getTangentSpace(centerMass, vertices, 2);

		tip = VectorDouble.linearCombination(1 + lambda, center, -lambda, centerMass);
		points = buildAngularPoints(vertices, angularSize, true);
		boundaryTangent = new VectorDouble[angularSize];
		centerTangent = new VectorDouble[angularSize];

		for (i = 0; i < angularSize; i++) {
			centerTangent[i] = VectorDouble.linearCombination(centerScale, points[i], -centerScale, centerMass);
			boundaryTangent[i] = VectorDouble.linearCombination(boundaryScale, points[i], -boundaryScale, tip);
			boundaryTangent[i] = VectorDouble.projection(tangentSpace, boundaryTangent[i], true);
		}

		nodes = new CellHigherNode[angularSize][radialSize];
		directed = new DefaultDirectedGraph<CellHigherNode, DefaultEdge>(DefaultEdge.class);

		for (j = 0; j < radialSize; j++) {
			t = ((double) (j + 1)) / ((double) radialSize);
			for (i = 0; i < angularSize; i++) {
				point = spline(t, center, points[i], centerTangent[i], boundaryTangent[i]);
				nodes[i][j] = new CellHigherNode(point, (j == radialSize - 1));
				directed.addVertex(nodes[i][j]);
			}
		}

		for (i = 0; i < angularSize; i++) {
			for (j = 0; j < radialSize; j++) {
				directed.addEdge(nodes[i][j], nodes[(i + 1) % (angularSize)][j]);
				if (j < radialSize - 1) {
					directed.addEdge(nodes[i][j], nodes[i][j + 1]);
				}
			}
		}

		if (placeCenter) {
			centerNode = new CellHigherNode(new VectorDouble(center), false);
			directed.addVertex(centerNode);

			for (i = 0; i < angularSize; i++) {
				directed.addEdge(centerNode, nodes[i][0]);
			}
		}

		return new AsUndirectedGraph<CellHigherNode, DefaultEdge>(directed);
	}
	
	/**
	 * Performs cubic spline interpolation.
	 * 
	 * @param t		how far to move from the first endpoint towards the second
	 * @param a		first endpoint
	 * @param b		second endpoint
	 * @param da	first velocity
	 * @param db	second velocity
	 * @return		the point
	 */
	public static VectorDouble spline(double t, VectorDouble a, VectorDouble b, VectorDouble da, VectorDouble db) {
		VectorDouble v2, v3, point;
		double t2, t3;
	
		t2 = t * t;
		t3 = t2 * t;
	
		v2 = VectorDouble.linearCombination(-3, a, 3, b, -2, da, -1, db);
		v3 = VectorDouble.linearCombination(2, a, -2, b, 1, da, 1, db);
		point = VectorDouble.linearCombination(1, a, t, da, t2, v2, t3, v3);
	
		return point;
	}

}
