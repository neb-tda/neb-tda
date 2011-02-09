package edu.stanford.math.nebtda;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.StringTokenizer;

import org.jgrapht.graph.DefaultEdge;

/**
 * Reads input and writes output.
 * 
 * @author Henry Adams
 * @author Atanas Atanasov
 */	
public class IO {
	
	/**
	 * Reads vector data from a file.
	 * 
	 * @param filename		name of the input file
	 * @return 				vectors in the file
	 */
	public static VectorDouble[] readDataFromFile(String filename) throws IOException {
		FileInputStream fis = new FileInputStream(filename);
		InputStreamReader isr = new InputStreamReader(fis);
		BufferedReader br = new BufferedReader(isr);

		ArrayList<ArrayList<Double>> datalist = new ArrayList<ArrayList<Double>>();
		ArrayList<Double> dataline;
		StringTokenizer tokenizer;
		String line;

		while ((line = br.readLine()) != null) {
			if (!line.trim().isEmpty()) {
				dataline = new ArrayList<Double>();
				datalist.add(dataline);
				tokenizer = new StringTokenizer(line);

				while (tokenizer.hasMoreTokens()) {
					dataline.add(Double.valueOf(tokenizer.nextToken().trim()));
				}
				tokenizer = null;
			}
		}
		if (datalist.size() == 0) {
			throw new IllegalArgumentException("No input lines read.");
		}

		br.close();
		isr.close();
		fis.close();

		int numPoints;
		numPoints = datalist.size();
		VectorDouble[] data = new VectorDouble[numPoints];

		for (int i = 0; i < numPoints; i++) {
			dataline = datalist.get(i);
			data[i] = new VectorDouble(dataline.size());
			for (int j = 0; j < dataline.size(); j++) {
				data[i].set(j, dataline.get(j).doubleValue());
			}
		}

		return data;
	}
	
	/**
	 * Writes 0-cells to a file.
	 * 
	 * @param cells			0-cells
	 * @param filename		name of the file
	 */
	public static void writeZeroCellsToFile(Cell0[] cells, String filename) throws IOException {
		PrintWriter printWriter = new PrintWriter(new File(filename));
		for (int i = 0; i < cells.length; i++) {
			printWriter.println(cells[i].point.toString(false, false));
		}
		printWriter.close();
	}
	
	/**
	 * Reads 0-cells from a file.
	 * 
	 * @param filename		name of the file
	 * @return				0-cells
	 */
	public static Cell0[] readZeroCellsFromFile(String filename) throws IOException {
		VectorDouble[] vectors = IO.readDataFromFile(filename);
		Cell0[] result = new Cell0[vectors.length];
		for (int i = 0; i < vectors.length; i++) {
			result[i] = new Cell0(vectors[i]);
		}
		return result;
	}
	
	/**
	 * Writes 1-cells to a file.
	 * 
	 * @param cells			1-cells
	 * @param filename		name of the file
	 */
	public static void writeOneCellsToFile(Cell1[] cells, String filename) throws IOException {
		PrintWriter printWriter = new PrintWriter(new File(filename));
		for (int i = 0; i < cells.length; i++) {
			for (int j = 0; j < cells[i].points.length; j++) {
				printWriter.println(cells[i].points[j].toString(false, false));
			}
		}
		printWriter.close();
	}
	
	/**
	 * Reads 1-cells from a file. The number of points in each band must be constant.
	 * 
	 * @param filename			name of the file
	 * @param numBandPoints		number of points in each band
	 * @return					1-cells
	 */
	public static Cell1[] readOneCellsFromFile(String filename, int numBandPoints) throws IOException {
		int i, j;
		
		VectorDouble[] points = IO.readDataFromFile(filename);
		if (points.length % numBandPoints != 0) {
			throw new IllegalArgumentException("The number of points in the file must be a multiple of the number of points in each band.");
		}
		int numBands = (int) (points.length / numBandPoints);
		
		VectorDouble[][] bands = new VectorDouble[numBands][numBandPoints];
		for (i = 0; i < numBands; i++) {
			bands[i] = new VectorDouble[numBandPoints];
			for (j = 0; j < numBandPoints; j++) {
				bands[i][j] = points[i*numBandPoints + j];
			}
		}
		
		Cell1[] result = new Cell1[bands.length];
		for (i = 0; i < bands.length; i++) {
			result[i] = new Cell1(bands[i]);
		}
		return result;
	}
	
	/**
	 * Writes a higher cell to a points file and a graph file.
	 * 
	 * @param cell			higher dimensional cells
	 * @param name			prefix of the file names
	 */
	public static void writeHigherCellToFile(CellHigher cell, String name) throws IOException {
		writeHigherCellPointsToFile(cell, name + "_" + cell.cellDimension + "points.txt");
		writeHigherCellGraphToFile(cell, name + "_" + cell.cellDimension + "graph.txt");
	}
	
	/**
	 * Writes a higher cell to a points file.
	 * 
	 * @param cell			higher cell
	 * @param filename		name of the file
	 */
	public static void writeHigherCellPointsToFile(CellHigher cell, String filename) throws IOException {
		PrintWriter printWriter = new PrintWriter(new File(filename));
		for (CellHigherNode node : cell.graph.vertexSet()) {
			printWriter.println(node.point.toString(false, false));
		}
		printWriter.close();
	}
	
	/**
	 * Writes a higher cell to a graph file.
	 * 
	 * @param cell			higher cell
	 * @param filename		name of the file
	 */
	public static void writeHigherCellGraphToFile(CellHigher cell, String filename) throws IOException {
		PrintWriter printWriter = new PrintWriter(new File(filename));
		for (DefaultEdge edge : cell.graph.edgeSet()) {
			printWriter.println(cell.graph.getEdgeSource(edge) + " " + cell.graph.getEdgeTarget(edge));
		}
		printWriter.close();
	}

}
