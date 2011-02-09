import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

import org.jgrapht.UndirectedGraph;
import org.jgrapht.graph.DefaultEdge;

import edu.stanford.math.nebtda.BoundaryLoop;
import edu.stanford.math.nebtda.Cell1;
import edu.stanford.math.nebtda.CellHigher;
import edu.stanford.math.nebtda.CellHigherNode;
import edu.stanford.math.nebtda.Density;
import edu.stanford.math.nebtda.DistributionNormal;
import edu.stanford.math.nebtda.IO;
import edu.stanford.math.nebtda.InitialBand;
import edu.stanford.math.nebtda.InitialBandEuclidean;
import edu.stanford.math.nebtda.InitialBandSphere;
import edu.stanford.math.nebtda.Parameters0;
import edu.stanford.math.nebtda.Parameters1;
import edu.stanford.math.nebtda.ParametersHigher;
import edu.stanford.math.nebtda.TrialHelper;
import edu.stanford.math.nebtda.VectorDouble;
import edu.stanford.math.nebtda.VectorDoubleComparator;

/**
 * Main class.
 * 
 * @author Henry Adams
 * @author Atanas Atanasov
 */
public class Main {
	
	/**
	 * Enum DataSet stores the standard deviation for each data set, as well as a boolean corresponding to 
	 * whether or not the data set is normalized to lie on a unit sphere.
	 */
	public enum DataSet {
		GeneExpressions		(1.35, false),
		RangeImagePatches	(0.35, true),
		OpticalFlowPatches	(0.30, true),
		OpticalImagePatches (0.20, true),
		SocialNetwork		(0.45, false);
		
		protected final double std;
		protected final boolean normalizedToSphere;
		/**
		 * Constructor.
		 * 
		 * @param std					standard deviation for the data set
		 * @param normalizedToSphere	true if the data set is normalized to lie on a unit sphere
		 */
		DataSet(double std, boolean normalizedToSphere) {
			this.std = std;
			this.normalizedToSphere = normalizedToSphere;
		}
	}
	
	/**
	 * Main.
	 * 
	 * @param args
	 */
	
	public static void main(String[] args) throws IOException {
		
		DataSet dataSet = DataSet.GeneExpressions;
		// The "SocialNetwork" data set is not included with this code.
		// If you are interested in obtaining this data, please email Henry at henrya@math.stanford.edu for further instructions.
		
		Density density = new Density(DistributionNormal.instance, dataSet.std, "data/" + dataSet + ".txt");
		TrialHelper t = new TrialHelper();
		
		//zeroCellTrials(dataSet, density, t);
		//oneCellTrials(dataSet, density, t);
		//twoCellTrial(dataSet, density, t);
		
		//zeroCellTrial(dataSet, density);
		oneCellTrial(dataSet, density, t);
		
	}
	
	/**
	 * Returns the parameters for 0-cell trials.
	 * 
	 * @param dataSet		data set
	 * @param density		density estimator
	 * @return				parameters for 0-cell trials
	 */
	public static Parameters0 getParameters0(DataSet dataSet, Density density) {
		double convergence = 0.0001;
		int numTrials = 30;
		double clusteringThreshold = 0.3;
		return new Parameters0(convergence, numTrials, clusteringThreshold);
	}
	
	/**
	 * Returns the parameters for 1-cell trials.
	 * 
	 * @param dataSet		data set
	 * @param density		density estimator
	 * @return				parameters for 1-cell trials
	 */
	public static Parameters1 getParameters1(DataSet dataSet, Density density) {
		double gradConst = 1.0 / density.maxGradNorm();
		double convergence = 0.0001;
		int numBandPoints = 11;
		InitialBand initialBand = dataSet.normalizedToSphere ? new InitialBandSphere(numBandPoints) : new InitialBandEuclidean(numBandPoints);
		double anglePnFnZero = Math.PI/6;
		double anglePnFnOne = Math.PI/2;
		double step = 0.1;
		int maxSteps = 6000;
		int numTrials = 2;
		int numSuccesses = 2;
		double nearnessThreshold = 0.5;
		double clusteringThreshold = 0.3;
		return new Parameters1(gradConst, convergence, initialBand, anglePnFnZero, anglePnFnOne, step, maxSteps, numTrials, numSuccesses, nearnessThreshold, clusteringThreshold);
	}
	
	/**
	 * Returns the parameters for 2-cell trials.
	 * 
	 * @param dataSet		data set
	 * @param density		density estimator
	 * @return				parameters for 0-cell trials
	 */
	public static ParametersHigher getParameters2(DataSet dataSet, Density density) {
		double gradConst = 1.0 / density.maxGradNorm();
		double convergence = 0.001;
		double step = 0.05;
		int maxSteps = 6000;
		return new ParametersHigher(gradConst, convergence, step, maxSteps);
	}
	
	/**
	 * Runs a single 0-cell trial using the mean shift method.
	 * 
	 * @param dataSet		data set
	 * @param density		density estimator
	 */
	public static void zeroCellTrial(DataSet dataSet, Density density) throws IOException {
		Parameters0 p0 = getParameters0(dataSet, density);
		TrialHelper.findZeroCell(density, p0, "data/output/" + dataSet);
	}
	
	/**
	 * Runs 0-cell trials using the mean shift method.
	 * 
	 * @param dataSet		data set
	 * @param density		density estimator
	 * @param t				trial helper
	 */
	public static void zeroCellTrials(DataSet dataSet, Density density, TrialHelper t) throws IOException {
		Parameters0 p0 = getParameters0(dataSet, density);
		t.findZeroCells(density, p0);
		
		switch(dataSet) {
		case OpticalImagePatches:
			t.thresholdZeroCells(0.08);
			break;
		case SocialNetwork:
			t.thresholdZeroCells(0.0006);
			break;
		}
		
		IO.writeZeroCellsToFile(t.getZeroCells(), "data/output/" + dataSet + "_0cells.txt");
		t.printZeroCellDensities();
		//t.printZeroCells();
		//t.printZeroCellDistances();
	}
	
	/**
	 * Runs a single 1-cell trial.
	 * 
	 * @param dataSet		data set
	 * @param density		density estimator
	 * @param t				trial helper
	 */
	public static void oneCellTrial(DataSet dataSet, Density density, TrialHelper t) throws IOException {
		Parameters1 p1 = getParameters1(dataSet, density);
		t.setZeroCells(IO.readZeroCellsFromFile("data/output/" + dataSet + "_0cells_ordered.txt"), density);
		TrialHelper.findOneCell(density, p1, t.getZeroCells()[0].getPoint(), t.getZeroCells()[1].getPoint(), "data/output/" + dataSet);
	}
	
	/**
	 * Runs 1-cell trials.
	 * 
	 * @param dataSet		data set
	 * @param density		density estimator
	 * @param t				trial helper
	 */
	public static void oneCellTrials(DataSet dataSet, Density density, TrialHelper t) throws IOException {		
		switch(dataSet) {
		case OpticalImagePatches:
			OpticalImagePatchesOneCellTrials(density, t);
			break;
		default:
			Parameters1 p1 = getParameters1(dataSet, density);
			t.setZeroCells(IO.readZeroCellsFromFile("data/output/" + dataSet + "_0cells_ordered.txt"), density);
			
			t.findOneCells(density, p1);
			break;
		}
		IO.writeOneCellsToFile(t.getOneCells(), "data/output/" + dataSet + "_1cells.txt");
		t.printOneCellDensities();
		//t.printOneCells();
		//t.printOneCellDistances();
	}
	
	/**
	 * Runs the 1-cell trials for the OpticalImagePatches data set. The primary circle 1-cells are easy to find,
	 * and are found in the same fashion as 1-cells in the other data sets. The secondary circle 1-cells are harder
	 * to find. To do so, we consider numAntipodalTrials initial bands between DCT basis vectors e1 and -e1 (resp. e2 and -e2).
	 * We sort these initial bands (roughly) based on which ones are most likely to flow to the secondary circles.
	 * We only run the first numAntipodalTrialsRun initial bands, generally finding the secondary circles. Of course, one could
	 * run all numAntipodalTrials bands, but this requires a lot of simulation time.
	 * 
	 * @param density		density estimator
	 * @param t				trial helper
	 */
	public static void OpticalImagePatchesOneCellTrials(Density density, TrialHelper t) throws IOException {
		Parameters1 p1 = getParameters1(DataSet.OpticalImagePatches, density);	
		t.setZeroCells(IO.readZeroCellsFromFile("data/output/" + DataSet.OpticalImagePatches + "_0cells_ordered.txt"), density);
		
		System.out.println("Primary circle trials:");
		for (int i = 0; i < 4; i++) {
			t.findOneCells(density, p1, i, (i + 1) % (4));
		}
		
		System.out.println("Secondary circle trials:");
		int numSecondaryTrials = 600;
		int numSecondaryTrialsRun = 6;
		for (int i = 0; i < 2; i++) {
			ArrayList<VectorDouble> dirs = new ArrayList<VectorDouble>(numSecondaryTrials);
			for (int j = 0; j < numSecondaryTrials; j++) {
				dirs.add(DistributionNormal.nextSphereVector(density.getDimension()));
			}	
			VectorDoubleComparator vdc = new VectorDoubleComparator(i + 2);
			Collections.sort(dirs, vdc);
			
			Cell1[] cells = new Cell1[numSecondaryTrialsRun];
			for (int j = 0; j < numSecondaryTrialsRun; j++) {
				InitialBandSphere initialBand = (InitialBandSphere) p1.getInitialBand();
				cells[j] = new Cell1(initialBand.getInitialBand(t.getZeroCells()[i].getPoint(), t.getZeroCells()[i + 2].getPoint(), dirs.get(j)));
			}
			
			t.findOneCells(density, p1, cells);
		}
	}

	/**
	 * Runs a single 2-cell trial.
	 * 
	 * @param dataSet		data set
	 * @param density		density estimator
	 * @param t				trial helper
	 */
	public static void twoCellTrial(DataSet dataSet, Density density, TrialHelper t) throws IOException {
		ParametersHigher p = getParameters2(dataSet, density);
		t.setOneCells(IO.readOneCellsFromFile("data/output/" + dataSet + "_1cells_ordered.txt", 11), density);

		Cell1[] boundaryCells;
		boolean[] boundaryOrientations;
		switch(dataSet) {
		case RangeImagePatches:
		case OpticalFlowPatches:
		case SocialNetwork:
			boundaryCells = new Cell1[] {t.getOneCells()[0], t.getOneCells()[2], t.getOneCells()[3], t.getOneCells()[1]};
			boundaryOrientations = new boolean[] {true, true, true, false};
			break;
		case GeneExpressions:
			boundaryCells = new Cell1[] {t.getOneCells()[0], t.getOneCells()[2], t.getOneCells()[1]};
			boundaryOrientations = new boolean[] {true, true, false};
			break;
		case OpticalImagePatches:
		default:	
			throw new IllegalArgumentException("We have not yet tested 2-cells in the optical image patches data set.");
			// TODO: test
		}
		BoundaryLoop loop = new BoundaryLoop(boundaryCells, boundaryOrientations);
		
		int meshDimension = 2;
		int radialSize = 10;
		int angularSize = 20;
		VectorDouble center = loop.computeAvgCenter();
		UndirectedGraph<CellHigherNode, DefaultEdge> graph = CellHigher.buildWeb(loop, angularSize, radialSize, center, true);		
		CellHigher cell = new CellHigher(graph, meshDimension, null);
		
		cell = TrialHelper.findHigherCell(density, p, cell, "data/output/" + dataSet);
		IO.writeHigherCellToFile(cell, "data/output/" + dataSet);
	}
	
}