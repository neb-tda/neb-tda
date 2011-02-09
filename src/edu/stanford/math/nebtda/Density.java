package edu.stanford.math.nebtda;

import java.io.FileNotFoundException;
import java.io.IOException;

/**
 * Implements a density estimator for data.
 * 
 * @author Henry Adams
 * @author Atanas Atanasov
 */
public class Density extends Potential{

	protected Distribution distribution;
	protected double std;
	protected VectorDouble[] data;
	protected int dimension;
	protected int numPoints;

	/**
	 * Constructor which accepts an array of VectorDoubles.
	 * 
	 * @param distribution	kernel used to estimate density
	 * @param std			standard deviation of the kernel
	 * @param data			data
	 */
	public Density(Distribution distribution, double std, VectorDouble[] data) {
		setDistribution(distribution);
		setStd(std);
		setData(data);
	}

	/**
	 * Constructor which reads the data from a file.
	 * 
	 * @param distribution	kernel used to estimate density
	 * @param std			standard deviation of the kernel
	 * @param filename		name of the input file
	 */
	public Density(Distribution distribution, double std, String filename) throws FileNotFoundException, IOException {
		setDistribution(distribution);
		setStd(std);
		setData(IO.readDataFromFile(filename));
	}

	/**
	 * Sets the distribution.
	 * 
	 * @param distribution	kernel used to estimate density
	 */
	public void setDistribution(Distribution distribution) {
		if (distribution == null) {
			throw new IllegalArgumentException("Distribution reference is null.");
		}
		this.distribution = distribution;
	}

	/**
	 * Gets the distribution.
	 * 
	 * @return distribution		kernel used to estimate density
	 */
	public Distribution getDistribution() {
		return distribution;
	}

	/**
	 * Sets the standard deviation.
	 * 
	 * @param std	standard deviation of the kernel
	 */
	public void setStd(double std) {
		if (std <= 0) {
			throw new IllegalArgumentException("Standard deviation must be positive.");
		}
		this.std = std;
	}

	/**
	 * Gets the standard deviation.
	 * 
	 * @return standard deviation of the kernel
	 */
	public double getStd() {
		return std;
	}

	/**
	 * Sets the data, the dimension, and the number of points.
	 * 
	 * @param data		data
	 */
	public void setData(VectorDouble[] data) {
		dimension = VectorDouble.checkConsistentDimension(data);
		numPoints = data.length;
		this.data = data;
	}

	/**
	 * Gets the data.
	 * 
	 * @return data
	 */
	public VectorDouble[] getData() {
		return data;
	}

	/**
	 * Gets the dimension of the data.
	 * 
	 * @return dimension of the data
	 */
	public int getDimension() {
		return dimension;
	}

	/**
	 * Gets the number of points in the data.
	 * 
	 * @return number of points in the data
	 */
	public int getNumPoints() {
		return numPoints;
	}

	/**
	 * Computes the density value at a point.
	 * 
	 * @param point		point
	 * @return 			density value
	 */
	public double getValue(VectorDouble point) {
		double result = 0;
		for (int i = 0; i < numPoints; i++) {
			result += distribution.getValue(point, data[i], std);
		}
		return result / numPoints;
	}

	/**
	 * Computes the gradient vector at a point.
	 * 
	 * @param point		point
	 * @return 			gradient vector
	 */
	public VectorDouble getGradient(VectorDouble point) {
		VectorDouble result = new VectorDouble(dimension);
		for (int i = 0; i < numPoints; i++) {
			result.add(distribution.getGradient(point, data[i], std));
		}
		result.scalarMult(1.0 / numPoints);
		return result;
	}
	
	/**
	 * Computes the density value and gradient vector at a point.
	 * 
	 * @param point		point
	 * @return 			density value and gradient vector
	 */
	public Tuple<Double, VectorDouble> getValueGradient(VectorDouble point) {
		VectorDouble gradient = new VectorDouble(dimension);
		Double value = new Double(0.0);
		Tuple<Double, VectorDouble> result;

		for (int i = 0; i < numPoints; i++) {
			result = distribution.getValueGradient(point, data[i], std);
			value += result.get1();
			gradient.add(result.get2());
		}

		value = value / numPoints;
		gradient.scalarMult(1.0 / numPoints);

		return new Tuple<Double, VectorDouble>(value, gradient);
	}
	
	/**
	 * Computes the maximum norm of all gradient vectors.
	 * 
	 * @return 			maximum norm of all gradient vectors
	 */
	public double maxGradNorm() {
		return distribution.maxGradNorm(std, dimension);
	}
	
	/**
	 * Computes the mean shift at a point. See equation (8) of "Mean shift, mode seeking, and clustering" 
	 * by Yizong Cheng.
	 * 
	 * @param point		point
	 * @return 			mean shift
	 */
	public VectorDouble getMeanShift(VectorDouble point) {
		VectorDouble result = new VectorDouble(dimension);
		double normalization = 0.0;
		double value;
		for (int i = 0; i < numPoints; i++) {
			value = distribution.getValue(point, data[i], std);
			result.linearCombination(1, value, data[i]);
			normalization += value;
		}
		result.scalarMult(1.0 / normalization);
		return result;
	}
	
}
