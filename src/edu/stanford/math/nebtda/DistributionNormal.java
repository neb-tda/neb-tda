package edu.stanford.math.nebtda;

import java.util.Random;

/**
 * Implements a normal distribution.
 * 
 * @author Henry Adams
 * @author Atanas Atanasov
 */
public class DistributionNormal extends Distribution {

	protected static Random r = new Random();
	public static DistributionNormal instance;
	static {
		instance = new DistributionNormal();
	}
	
	/**
	 * Empty constructor.
	 */
	public DistributionNormal() {
	}

	/**
	 * Computes the value of the distribution at a point.
	 * 
	 * @param point		point
	 * @param mean		mean of the distribution
	 * @param std		standard deviation of the distribution
	 * @return			value at the point
	 */
	@Override
	public double getValue(VectorDouble point, VectorDouble mean, double std) {
		return Math.exp(-VectorDouble.subtract(mean, point).normSquared() / (2 * std * std)) / Math.pow(Math.sqrt(2 * Math.PI) * std, point.coordinates.length);
	}

	/**
	 * Computes the gradient of the distribution at a point.
	 * 
	 * @param point		point
	 * @param mean		mean of the distribution
	 * @param std		standard deviation of the distribution
	 * @return			gradient at the point
	 */
	@Override
	public VectorDouble getGradient(VectorDouble point, VectorDouble mean, double std) {
		VectorDouble result = VectorDouble.subtract(mean, point);
		result.scalarMult(Math.exp(-result.normSquared() / (2 * std * std)) / (Math.pow(Math.sqrt(2 * Math.PI) * std, point.coordinates.length) * std * std));
		return result;
	}

	/**
	 * Computes the value and gradient of the distribution at a point.
	 * 
	 * @param point		point
	 * @param mean		mean of the distribution
	 * @param std		standard deviation of the distribution
	 * @return			value and gradient at the point
	 */
	@Override
	public Tuple<Double, VectorDouble> getValueGradient(VectorDouble point, VectorDouble mean, double std) {
		VectorDouble result = VectorDouble.subtract(mean, point);
		double value = Math.exp(-result.normSquared() / (2 * std * std)) / Math.pow(Math.sqrt(2 * Math.PI) * std, point.coordinates.length);
		result.scalarMult(value / (std * std));
		return new Tuple<Double, VectorDouble>(new Double(value), result);
	}
	
	/**
	 * Computes the maximum norm of the gradient of the distribution.
	 * 
	 * @param std			standard deviation of the distribution
	 * @param dimension		dimension of the distribution
	 * @return				maximum norm of the gradient
	 */
	public double maxGradNorm(double std, int dimension) {
		return 1.0 / (Math.pow(Math.sqrt(2 * Math.PI) * std, dimension) * Math.sqrt(Math.E));
	}
	
	/**
	 * Returns a vector from the normal distribution.
	 * 
	 * @param dimension		dimension of the vector
	 * @return				a normally distributed vector
	 */
	public static VectorDouble nextNormalVector(int dimension) {
		if (dimension < 1) {
			throw new IllegalArgumentException("Dimension must be positive.");
		}
		
		double[] normals = new double[dimension];
		for (int i = 0; i < dimension; i++) {
			normals[i] = r.nextGaussian();
		}
		
		return new VectorDouble(normals);
	}
	
	/**
	 * Returns a vector from the uniform distribution on the sphere. The specified dimension is the dimension
	 * of the ambient Euclidean space, which is one more than the dimension of the sphere as a manifold.
	 * 
	 * @param dimension		dimension
	 * @return				random sphere point
	 */
	public static VectorDouble nextSphereVector(int dimension) {
		if (dimension < 1) {
			throw new IllegalArgumentException("Dimension must be positive.");
		}
		
		VectorDouble point;
		double norm;
	
		for (;;) {
			point = DistributionNormal.nextNormalVector(dimension);
			norm = point.norm();
			if (norm > 0) {
				point.scalarMult(1.0 / norm);
				return point;
			}
		}
	}

}