package edu.stanford.math.nebtda;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Collection;
import java.util.Locale;
import java.util.Random;

import org.jgrapht.graph.DefaultEdge;

import edu.stanford.math.nebtda.Cell0;
import edu.stanford.math.nebtda.Cell1;
import edu.stanford.math.nebtda.CellHigher;
import edu.stanford.math.nebtda.Density;
import edu.stanford.math.nebtda.Parameters0;
import edu.stanford.math.nebtda.VectorDouble;

/**
 * Contains the tools for running trials.
 * 
 * @author Henry Adams
 * @author Atanas Atanasov
 */
public class TrialHelper {
	
	protected static PrintWriter nullPrintWriter = new PrintWriter(NullOutputStream.NULL_OUTPUT_STREAM);
	protected static Random r = new Random();
	protected static int SMALL_DATA_SET_SIZE = 100;

	protected Cell0[] zeroCells;
	protected Cell1[] oneCells;
	protected CellHigher[] higherCells;
	
	protected DecimalFormat decSci = new DecimalFormat("0.###E0;-0.###E0", new DecimalFormatSymbols(Locale.US));
	
	/**
	 * Constructor.
	 */
	public TrialHelper() {
		zeroCells = new Cell0[0];
		oneCells = new Cell1[0];
		higherCells = new CellHigher[0];
	}
	
	/**
	 * Sets the 0-cells.
	 * 
	 * @param zeroCells		0-cells
	 */
	public void setZeroCells(Cell0[] zeroCells) {
		this.zeroCells = zeroCells;
	}
	
	/**
	 * Sets the 1-cells.
	 * 
	 * @param oneCells		1-cells
	 */
	public void setOneCells(Cell1[] oneCells) {
		this.oneCells = oneCells;
	}
	
	/**
	 * Sets the higher cells.
	 * 
	 * @param higherCells		higher cells
	 */
	public void setHigherCells(CellHigher[] higherCells) {
		this.higherCells = higherCells;
	}
	
	/**
	 * Sets the 0-cells and their density values.
	 * 
	 * @param zeroCells		0-cells
	 * @param density		density estimator
	 */
	public void setZeroCells(Cell0[] zeroCells, Density density) {
		computeCellDensities(zeroCells, density);
		setZeroCells(zeroCells);
	}
	
	/**
	 * Sets the 1-cells and their density values.
	 * 
	 * @param oneCells		1-cells
	 * @param density		density estimator
	 */
	public void setOneCells(Cell1[] oneCells, Density density) {
		computeCellDensities(oneCells, density);
		setOneCells(oneCells);
	}
	
	/**
	 * Sets the higher cells and their density values.
	 * 
	 * @param higherCells		higherCells
	 * @param density			density estimator
	 */
	public void setHigherCells(CellHigher[] higherCells, Density density) {
		computeCellDensities(higherCells, density);
		setHigherCells(higherCells);
	}
	
	/**
	 * Sets the density values of the cells.
	 * 
	 * @param cells			cells
	 * @param density		density estimator
	 */
	public void computeCellDensities(Cell[] cells, Density density) {
		for (int i = 0; i < cells.length; i++) {
			cells[i].setValue(cells[i].computeValue(density));
		}
	}
	
	/**
	 * Gets the 0-cells.
	 * 
	 * @return 0-cells
	 */
	public Cell0[] getZeroCells() {
		return zeroCells;
	}
	
	/**
	 * Gets the 1-cells.
	 * 
	 * @return 1-cells
	 */
	public Cell1[] getOneCells() {
		return oneCells;
	}
	
	/**
	 * Gets the higher cells.
	 * 
	 * @return higher cells
	 */
	public CellHigher[] getHigherCells() {
		return higherCells;
	}
	
	/**
	 * Adds new 0-cells.
	 * 
	 * @param newCells		0-cells to add
	 */
	public void addZeroCells(Cell0[] newCells) {
		int i;
		Cell0[] cells = new Cell0[zeroCells.length + newCells.length];
		for (i = 0; i < zeroCells.length; i++) {
			cells[i] = zeroCells[i];
		}
		for (i = 0; i < newCells.length; i++) {
			cells[zeroCells.length + i] = newCells[i];
		}
		this.zeroCells = cells;
	}
	
	/**
	 * Adds new 1-cells.
	 * 
	 * @param newCells		1-cells to add
	 */
	public void addOneCells(Cell1[] newCells) {
		int i;
		Cell1[] cells = new Cell1[oneCells.length + newCells.length];
		for (i = 0; i < oneCells.length; i++) {
			cells[i] = oneCells[i];
		}
		for (i = 0; i < newCells.length; i++) {
			cells[oneCells.length + i] = newCells[i];
		}
		this.oneCells = cells;
	}
	
	/**
	 * Adds new higher cells.
	 * 
	 * @param newCells		higher cells to add
	 */
	public void addHigherCells(CellHigher[] newCells) {
		int i;
		CellHigher[] cells = new CellHigher[higherCells.length + newCells.length];
		for (i = 0; i < higherCells.length; i++) {
			cells[i] = higherCells[i];
		}
		for (i = 0; i < newCells.length; i++) {
			cells[higherCells.length + i] = newCells[i];
		}
		this.higherCells = cells;
	}
	
	/**
	 * Searches for a single 0-cell using the mean shift method.
	 * 
	 * @param density 		density estimator
	 * @param p0			trial parameters
	 * @param name			name of the trial
	 * @return 				0-cell
	 */
	public static Cell0 findZeroCell(Density density, Parameters0 p0, String name) throws IOException {
		PrintWriter pointWriter = new PrintWriter(new File(name + "_0points.txt"));
		PrintWriter distanceWriter = new PrintWriter(new File(name + "_0forces.txt"));

		Cell0 cell = new Cell0(new VectorDouble(density.data[r.nextInt(density.getNumPoints())]));
		cell.run(density, p0, true, pointWriter, distanceWriter);
		System.out.println("density " + cell.value);
		
		pointWriter.close();
		distanceWriter.close();
		return cell;
	}
	
	/**
	 * Searches for a single 1-cell. If the 1-cell did not converge, then its points are null.
	 * 
	 * @param density 		density estimator
	 * @param p1			trial parameters
	 * @param start			first endpoint
	 * @param end			second endpoint
	 * @param name			name of the trial
	 * @return 				1-cell
	 */
	public static Cell1 findOneCell(Density density, Parameters1 p1, VectorDouble start, VectorDouble end, String name) throws IOException {	
		PrintWriter pointWriter = new PrintWriter(new File(name + "_1points.txt"));
		PrintWriter forceWriter = new PrintWriter(new File(name + "_1forces.txt"));
		
		Cell1 cell = new Cell1(p1.initialBand.getInitialBand(start, end));	
		cell.run(density, p1, true, pointWriter, forceWriter);
		if (cell.converged == false) {
			System.out.println("No convergence");
		} else {
			System.out.println("density " + cell.value);
		}
		
		pointWriter.close();
		forceWriter.close();
		return cell;
	}

	/**
	 * Searches for a single higher cell. If the higher cell did not converge, then its points are null.
	 * 
	 * @param density 		density estimator
	 * @param p				trial parameters
	 * @param cell			cell
	 * @param name			name of the trial
	 * @return 				higher cell cell
	 */
	public static CellHigher findHigherCell(Density density, ParametersHigher p, CellHigher cell, String name) throws IOException {
		PrintWriter pointWriter = nullPrintWriter;
		PrintWriter forceWriter = new PrintWriter(new File(name + "_" + cell.getMeshDimension() + "forces.txt"));
		
		cell.run(density, p, true, pointWriter, forceWriter);
		if (cell.converged == false) {
			System.out.println("No convergence");
		} else {
			System.out.println("density " + cell.value);
		}
		
		pointWriter.close();
		forceWriter.close();
		return cell;
	}
	
	/**
	 * Searches for 0-cells using the mean shift method.
	 * 
	 * @param density 		density estimator
	 * @param p0			trial parameters
	 */
	public void findZeroCells(Density density, Parameters0 p0) {
		Cell0[] vertices;
		int i;
		if (density.getNumPoints() > SMALL_DATA_SET_SIZE) {
			
			vertices = new Cell0[p0.numTrials];
			for (i = 0; i < p0.numTrials; i++) {
				System.out.println("Trial " + (i + 1) + "...");
				Cell0 cell = new Cell0(new VectorDouble(density.data[r.nextInt(density.getNumPoints())]));
				cell.run(density, p0, false, nullPrintWriter, nullPrintWriter);
				vertices[i] = cell;
			}			
			
		} else {
			
			vertices = new Cell0[density.getNumPoints()];
			for (i = 0; i < density.getNumPoints(); i++) {
				System.out.println("Initial data point " + (i + 1) + "...");
				Cell0 cell = new Cell0(new VectorDouble(density.data[i]));
				cell.run(density, p0, false, nullPrintWriter, nullPrintWriter);
				vertices[i] = cell;
			}
		}
		addZeroCells(clusterZeroCells(vertices, p0.clusteringThreshold));
	}
	
	/**
	 * Searches for 1-cells.
	 * 
	 * @param density 		density estimator
	 * @param p1			trial parameters
	 */
	public void findOneCells(Density density, Parameters1 p1) {
		for (int i = 0; i < zeroCells.length; i++) {
			for (int j = i + 1; j < zeroCells.length; j++) {
				findOneCells(density, p1, i, j);
			}
		}
	}
	
	/**
	 * Searches for 1-cells between two given endpoints.
	 * 
	 * @param density 		density estimator
	 * @param p1			trial parameters
	 * @param index1		index of the first 0-cell endpoint
	 * @param index2		index of the second 0-cell endpoint
	 */
	public void findOneCells(Density density, Parameters1 p1, int index1, int index2) {
		if ((index1 < 0) || (index2 < 0) || (index1 >= zeroCells.length) || (index2 >= zeroCells.length)) {
			throw new IllegalArgumentException("The indices must be between 0 and the number of zero cells minus one.");
		}
		if (index1 == index2) {
			throw new IllegalArgumentException("The indices must not be equal.");
		}
		
		Cell1 cell;
		int i;
		
		VectorDouble[] vertices = new VectorDouble[zeroCells.length];
		for (i = 0; i < vertices.length; i++) {
			vertices[i] = new VectorDouble(zeroCells[i].point);
		}
		
		System.out.println("Trials between 0-cells " + (index1 + 1) + " and " + (index2 + 1) + ":");
		
		VectorDouble[] otherVertices = new VectorDouble[vertices.length - 2];
		int vertexCount = 0;
		for (i = 0; i < vertices.length; i++) {
			if ((i != index1) && (i != index2)) {
				otherVertices[vertexCount] = vertices[i];
				vertexCount ++;
			}
		}
		
		LinkedList<Cell1> bands = new LinkedList<Cell1>();
		int successCount = 0;
		for (i = 0; (successCount < p1.numSuccesses) && (i < p1.numTrials); i++) {
			System.out.print("Trial " + (i + 1) + ": ");
			cell = new Cell1(p1.initialBand.getInitialBand(vertices[index1], vertices[index2]));
			cell.run(density, p1, false, nullPrintWriter, nullPrintWriter);
			if (cell.converged == false) {
				System.out.println("no convergence");
			} else if (cell.bandNearVertices(otherVertices, p1.nearnessThreshold)) {
				System.out.println("near other vertices");
			} else {
				bands.add(cell);
				successCount++;
				System.out.println("success number " + successCount);
			}
		}
		addOneCells(clusterOneCells(bands.toArray(new Cell1[0]), p1.clusteringThreshold));
	}
	
	/**
	 * Runs the given initial 1-cells.
	 * 
	 * @param density 		density estimator
	 * @param p1			trial parameters
	 * @param cells			initial 1-cells
	 */
	public void findOneCells(Density density, Parameters1 p1, Cell1[] cells) {
		Cell1 cell;
		LinkedList<Cell1> bands = new LinkedList<Cell1>();
		
		int successCount = 0;
		for(int i = 0; i < cells.length; i++) {
			System.out.print("Trial " + (i + 1) + ": ");
			cell = cells[i];
			cell.run(density, p1, false, nullPrintWriter, nullPrintWriter);
			if (cell.converged == false) {
				System.out.println("no convergence");
			} else {
				bands.add(cell);
				successCount++;
				System.out.println("success number " + successCount);
			}
		}
		addOneCells(clusterOneCells(bands.toArray(new Cell1[0]), p1.clusteringThreshold));
	}
	
	/**
	 * Clusters 0-cells, using single-linkage clustering. Returns the densest 0-cell from each cluster.
	 * 
	 * @param vertices 				0-cells to cluster
	 * @param clusteringThreshold	clustering threshold
	 */
	public Cell0[] clusterZeroCells(Cell0[] vertices, double clusteringThreshold) {
		ClusteringSingleLinkage<Cell0> slc = new ClusteringSingleLinkage<Cell0>(Cell0.instance, clusteringThreshold);
		Collection<Cluster<Cell0>> collection = slc.cluster(vertices);

		ClustRepMaxValue<Cell0> crmv = new ClustRepMaxValue(Cell0.instance);
		Cell0[] cells = new Cell0[collection.size()];
		int i = 0;
		for (Cluster<Cell0> cluster : collection) {
			cells[i++] = crmv.getRepresentative(cluster);
		}
		return cells;
	}
	
	/**
	 * Clusters 1-cells, using single-linkage clustering. Returns the densest 1-cell from each cluster.
	 * 
	 * @param bands 				1-cells to cluster
	 * @param clusteringThreshold	clustering threshold
	 */
	public Cell1[] clusterOneCells(Cell1[] bands, double clusteringThreshold) {
		ClusteringSingleLinkage<Cell1> slc = new ClusteringSingleLinkage<Cell1>(Cell1.instance, clusteringThreshold);
		Collection<Cluster<Cell1>> collection = slc.cluster(bands);

		ClustRepMaxValue<Cell1> crmv = new ClustRepMaxValue(Cell0.instance);
		Cell1[] cells = new Cell1[collection.size()];
		int i = 0;
		for (Cluster<Cell1> cluster : collection) {
			//System.out.println("cluster size: " + cluster.size());
			cells[i++] = crmv.getRepresentative(cluster);
		}
		return cells;
	}
	
	/**
	 * Removes the 0-cells with density below a given threshold.
	 * 
	 * @param densityThreshold		density threshold value
	 */
	public void thresholdZeroCells(double densityThreshold) {
		ArrayList<Cell0> cellList = new ArrayList<Cell0>();
		for (int i = 0; i < zeroCells.length; i++) {
			if(zeroCells[i].value >= densityThreshold) {
				cellList.add(zeroCells[i]);
			}
		}
		zeroCells = cellList.toArray(new Cell0[0]);
	}
	
	/**
	 * Removes the 1-cells with density below a given threshold.
	 * 
	 * @param densityThreshold		density threshold value
	 */
	public void thresholdOneCells(double densityThreshold) {
		ArrayList<Cell1> cellList = new ArrayList<Cell1>();
		for (int i = 0; i < oneCells.length; i++) {
			if(oneCells[i].value >= densityThreshold) {
				cellList.add(oneCells[i]);
			}
		}
		oneCells = cellList.toArray(new Cell1[0]);
	}
	
	/**
	 * Removes the higher cells with density below a given threshold.
	 * 
	 * @param densityThreshold		density threshold value
	 */
	public void thresholdHigherCells(double densityThreshold) {
		ArrayList<CellHigher> cellList = new ArrayList<CellHigher>();
		for (int i = 0; i < higherCells.length; i++) {
			if(higherCells[i].value >= densityThreshold) {
				cellList.add(higherCells[i]);
			}
		}
		higherCells = cellList.toArray(new CellHigher[0]);
	}
	
	/**
	 * Prints the 0-cells.
	 */
	public void printZeroCells() {
		System.out.println("0-cells:");
		for (int i = 0; i < zeroCells.length; i++) {
			System.out.println((i + 1) + ": " + "density " + decSci.format(zeroCells[i].value));
			System.out.println(zeroCells[i].point.toString(false, false));
		}
	}
	
	/**
	 * Prints the 1-cells.
	 */
	public void printOneCells() {
		System.out.println("1-cells:");
		for (int i = 0; i < oneCells.length; i++) {
			System.out.println((i + 1) + ": " + "density " + decSci.format(oneCells[i].value));
			for (int j = 0; j < oneCells[i].points.length; j++) {
				System.out.println(oneCells[i].points[j].toString(false, false));
			}
		}
	}
	
	/**
	 * Prints the higher cells.
	 */
	public void printHigherCells() {
		System.out.println("Higher-cells:");
		for (int i = 0; i < higherCells.length; i++) {
			System.out.println((i + 1) + ": " + "density " + decSci.format(higherCells[i].value));
			System.out.println("points:");
			for (CellHigherNode node : higherCells[i].graph.vertexSet()) {
				System.out.println(node.point.toString(false, false));
			}
			System.out.println("graph:");
			for (DefaultEdge edge : higherCells[i].graph.edgeSet()) {
				System.out.println(higherCells[i].graph.getEdgeSource(edge) + " " + higherCells[i].graph.getEdgeTarget(edge));
			}
		}
	}
	
	/**
	 * Prints the densities of the 0-cells.
	 */
	public void printZeroCellDensities() {
		System.out.println("0-cell densities:");
		for (int i = 0; i < zeroCells.length; i++) {
			System.out.println((i + 1) + ": " + "density " + decSci.format(zeroCells[i].value));
		}
	}
	
	/**
	 * Prints the densities of the 1-cells.
	 */
	public void printOneCellDensities() {
		System.out.println("1-cell densities:");
		for (int i = 0; i < oneCells.length; i++) {
			System.out.println((i + 1) + ": " + "density " + decSci.format(oneCells[i].value));
		}
	}
	
	/**
	 * Prints the densities of the higher cells.
	 */
	public void printHigherCellDensities() {
		System.out.println("Higher-cell densities:");
		for (int i = 0; i < higherCells.length; i++) {
			System.out.println((i + 1) + ": " + "density " + decSci.format(higherCells[i].value));
		}
	}
	
	/**
	 * Prints the distances between all of the 0-cells.
	 */
	public void printZeroCellDistances() {
		int i, j;
		System.out.println("O-cell distances:");
		for (i = 0; i < zeroCells.length; i++) {
			for (j = i+1; j < zeroCells.length; j++) {
				System.out.println((i + 1) + " and "  + (j + 1) + ": distance " + decSci.format(Cell0.instance.distance(zeroCells[i], zeroCells[j])));
			}
		}
	}
	
	/**
	 * Prints the distances between all of the 0-cells.
	 */
	public void printOneCellDistances() {
		int i, j;
		System.out.println("1-cell distances:");
		for (i = 0; i < oneCells.length; i++) {
			for (j = i+1; j < oneCells.length; j++) {
				System.out.println((i + 1) + " and "  + (j + 1) + ": distance " + decSci.format(Cell1.instance.distance(oneCells[i], oneCells[j])));
			}
		}
	}
		
}