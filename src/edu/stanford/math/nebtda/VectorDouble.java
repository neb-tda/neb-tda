package edu.stanford.math.nebtda;

import java.io.StringWriter;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;

import Jama.Matrix;
import Jama.SingularValueDecomposition;

/**
 * Implements vectors in Euclidean space.
 * 
 * @author Henry Adams
 * @author Atanas Atanasov
 */
public class VectorDouble {
	
	protected double[] coordinates;

	/**
	 * Constructor which accepts an array of doubles.
	 * 
	 * @param coordinates	coordinates of the vector	
	 */
	public VectorDouble(double[] coordinates) {
		setCoordinates(coordinates);
	}

	/**
	 * Constructor which copies another vector.
	 * 
	 * @param vector	vector to copy
	 */
	public VectorDouble(VectorDouble vector) {
		double[] oldCoordinates = vector.getCoordinates();
		double[] newCoordinates = new double[oldCoordinates.length];
		for (int i = 0; i < newCoordinates.length; i++) {
			newCoordinates[i] = oldCoordinates[i];
		}
		setCoordinates(newCoordinates);
	}

	/**
	 *  Constructor which creates the zero vector.
	 * 
	 * @param dimension		dimension of the zero vector
	 */
	public VectorDouble(int dimension) {
		this(new double[dimension]);
	}

	/**
	 * Sets the coordinates.
	 * 
	 * @param coordinates	coordinates
	 */
	public void setCoordinates(double[] coordinates) {
		if (coordinates == null) {
			throw new IllegalArgumentException("Coordinates reference is null.");
		}
		if (coordinates.length <= 0) {
			throw new IllegalArgumentException("Number of coordinates must be a positive integer.");
		}
		this.coordinates = coordinates;
	}
	
	/**
	 * Gets the coordinates.
	 * 
	 * @return coordinates
	 */
	public double[] getCoordinates() {
		return coordinates;
	}
	
	/**
	 * Sets the coordinate value at the index.
	 * 
	 * @param index			index
	 * @param value			value
	 */
	public void set(int index, double value) {
		if ((index >= 0) && (index < coordinates.length)) {
			coordinates[index] = value;
		} else {
			throw new IllegalArgumentException("Index is out of range.");
		}
	}

	/**
	 * Gets the coordinate value at the index.
	 * 
	 * @return coordinate
	 */
	public double get(int index) {
		if ((index >= 0) && (index < coordinates.length)) {
			return coordinates[index];
		} else {
			throw new IllegalArgumentException("Index is out of range.");
		}
	}

	public int getDimension() {
		return coordinates.length;
	}

	/**
	 * Assigns this vector to be the input vector.
	 * 
	 * @param a		input vector
	 */
	public void assign(VectorDouble a) {
		if (coordinates.length != a.coordinates.length) {
			throw new IllegalArgumentException("Input vectors do not have the same length.");
		}

		for (int i = 0; i < coordinates.length; i++) {
			coordinates[i] = a.coordinates[i];
		}
	}

	/**
	 * Checks if this vector is equal to another.
	 * 
	 * @param a		input vector
	 * @return		true if equal
	 */
	public boolean equals(VectorDouble a) {
		return equals(this, a);
	}

	/**
	 * Checks if two vectors are equal.
	 * 
	 * @param a		first input vector
	 * @param b		second input vector
	 * @return		true if equal
	 */
	public static boolean equals(VectorDouble a, VectorDouble b) {
		if (a.coordinates.length != b.coordinates.length) {
			throw new IllegalArgumentException("The input vectors do not have the same length.");
		}

		if (a == b) {
			return true;
		}
		for (int i = 0; i < a.coordinates.length; i++) {
			if (a.coordinates[i] != b.coordinates[i]) {
				return false;
			}
		}
		return true;
	}

	/**
	 * Adds a vector.
	 * 
	 * @param a		vector to add
	 */
	public void add(VectorDouble a) {
		if (this.coordinates.length != a.coordinates.length) {
			throw new IllegalArgumentException("Input vector does not have the same length.");
		}
		
		for (int i = 0; i < coordinates.length; i++) {
			coordinates[i] += a.coordinates[i];
		}
	}

	/**
	 * Adds two vectors.
	 * 
	 * @param a		first vector
	 * @param b		second vector
	 * @return		sum of the vectors
	 */
	public static VectorDouble add(VectorDouble a, VectorDouble b) {
		if (a.coordinates.length != b.coordinates.length) {
			throw new IllegalArgumentException("Input vectors do not have the same length.");
		}

		VectorDouble sum = new VectorDouble(a.coordinates.length);
		for (int i = 0; i < sum.coordinates.length; i++) {
			sum.coordinates[i] = a.coordinates[i] + b.coordinates[i];
		}
		return sum;
	}

	/**
	 * Subtracts a vector.
	 * 
	 * @param a		vector to subtract
	 */
	public void subtract(VectorDouble a) {
		if (this.coordinates.length != a.coordinates.length) {
			throw new IllegalArgumentException("The input vector does not have the same length.");
		}

		for (int i = 0; i < coordinates.length; i++) {
			coordinates[i] -= a.coordinates[i];
		}
	}

	/**
	 * Subtracts two vectors.
	 * 
	 * @param a		first vector
	 * @param b		second vector
	 * @return		vector a minus vector b
	 */
	public static VectorDouble subtract(VectorDouble a, VectorDouble b) {
		if (a.coordinates.length != b.coordinates.length) {
			throw new IllegalArgumentException("The input vectors do not have the same length.");
		}

		VectorDouble difference = new VectorDouble(a.coordinates.length);
		for (int i = 0; i < difference.coordinates.length; i++) {
			difference.coordinates[i] = a.coordinates[i] - b.coordinates[i];
		}
		return difference;
	}
	
	/**
	 * Multiplies this vector by a scalar.
	 * 
	 * @param c		scalar
	 */
	public void scalarMult(double c) {
		for (int i = 0; i < coordinates.length; i++) {
			coordinates[i] = c * coordinates[i];
		}
	}

	/**
	 * Multiplies a vector by a scalar.
	 * 
	 * @param c		scalar
	 * @param a		vector
	 * @return		scalar c times vector a
	 */
	public static VectorDouble scalarMult(double c, VectorDouble a) {
		VectorDouble product = new VectorDouble(a.coordinates.length);
		for (int i = 0; i < product.coordinates.length; i++) {
			product.coordinates[i] = c * a.coordinates[i];
		}
		return product;
	}

	/**
	 * Negates this vector.
	 */
	public void negate() {
		this.scalarMult(-1);
	}

	/**
	 * Negates a vector.
	 * 
	 * @param a		vector
	 * return		negated vector
	 */
	public static VectorDouble negate(VectorDouble a) {
		return scalarMult(-1, a);
	}

	/**
	 * Computes a linear combination.
	 * 
	 * @param c1	first scalar
	 * @param c2	second scalar
	 * @param a2	second vector
	 * return		(c1 times this) plus (c2 times a2)
	 */
	public void linearCombination(double c1, double c2, VectorDouble a2) {
		if (coordinates.length != a2.coordinates.length) {
			throw new IllegalArgumentException("Input vector does not have the same length.");
		}
		for (int i = 0; i < this.coordinates.length; i++) {
			coordinates[i] = c1 * coordinates[i] + c2 * a2.coordinates[i];
		}
	}

	/**
	 * Computes a linear combination.
	 * 
	 * @param c1	first scalar
	 * @param c2	second scalar
	 * @param a2	second vector
	 * @param c3	third scalar
	 * @param a3	third vector
	 * return		(c1 times this) plus (c2 times a2) plus (c3 times a3)
	 */
	public void linearCombination(double c1, double c2, VectorDouble a2, double c3, VectorDouble a3) {
		if ((coordinates.length != a2.coordinates.length) || (a2.coordinates.length != a3.coordinates.length)) {
			throw new IllegalArgumentException("Input vectors do not have the same length.");
		}
		for (int i = 0; i < this.coordinates.length; i++) {
			coordinates[i] = c1 * coordinates[i] + c2 * a2.coordinates[i] + c3 * a3.coordinates[i];
		}
	}

	/**
	 * Computes a linear combination.
	 * 
	 * @param c1	first scalar
	 * @param c2	second scalar
	 * @param a2	second vector
	 * @param c3	third scalar
	 * @param a3	third vector
	 * @param c4	fourth scalar
	 * @param a4	fourth vector
	 * return		(c1 times this) plus (c2 times a2) plus (c3 times a3) plus (c4 times a4)
	 */
	public void linearCombination(double c1, double c2, VectorDouble a2, double c3, VectorDouble a3, double c4,
			VectorDouble a4) {
		if ((coordinates.length != a2.coordinates.length) || (a2.coordinates.length != a3.coordinates.length)
				|| (a3.coordinates.length != a4.coordinates.length)) {
			throw new IllegalArgumentException("Input vectors do not have the same length.");
		}
		for (int i = 0; i < this.coordinates.length; i++) {
			coordinates[i] = c1 * coordinates[i] + c2 * a2.coordinates[i] + c3 * a3.coordinates[i] + c4
					* a4.coordinates[i];
		}
	}

	/**
	 * Computes a linear combination.
	 * 
	 * @param c1	first scalar
	 * @param a1	first vector
	 * @param c2	second scalar
	 * @param a2	second vector
	 * return		(c1 times a1) plus (c2 times a2)
	 */
	public static VectorDouble linearCombination(double c1, VectorDouble a1, double c2, VectorDouble a2) {
		if (a1.coordinates.length != a2.coordinates.length) {
			throw new IllegalArgumentException("Input vectors do not have the same length.");
		}
		VectorDouble linearCombo = new VectorDouble(a1.coordinates.length);
		for (int i = 0; i < linearCombo.coordinates.length; i++) {
			linearCombo.coordinates[i] = c1 * a1.coordinates[i] + c2 * a2.coordinates[i];
		}
		return linearCombo;
	}

	/**
	 * Computes a linear combination.
	 * 
	 * @param c1	first scalar
	 * @param a1	first vector
	 * @param c2	second scalar
	 * @param a2	second vector
	 * @param c3	third scalar
	 * @param a3	third vector
	 * return		(c1 times a1) plus (c2 times a2) plus (c3 times a3)
	 */
	public static VectorDouble linearCombination(double c1, VectorDouble a1, double c2, VectorDouble a2, double c3,
			VectorDouble a3) {
		if ((a1.coordinates.length != a2.coordinates.length) || (a2.coordinates.length != a3.coordinates.length)) {
			throw new IllegalArgumentException("Input vectors do not have the same length.");
		}
		VectorDouble linearCombo = new VectorDouble(a1.coordinates.length);
		for (int i = 0; i < linearCombo.coordinates.length; i++) {
			linearCombo.coordinates[i] = c1 * a1.coordinates[i] + c2 * a2.coordinates[i] + c3 * a3.coordinates[i];
		}
		return linearCombo;
	}

	/**
	 * Computes a linear combination.
	 * 
	 * @param c1	first scalar
	 * @param a1	first vector
	 * @param c2	second scalar
	 * @param a2	second vector
	 * @param c3	third scalar
	 * @param a3	third vector
	 * @param c4	fourth scalar
	 * @param a4	fourth vector
	 * return		(c1 times a1) plus (c2 times a2) plus (c3 times a3) plus (c4 times a4)
	 */
	public static VectorDouble linearCombination(double c1, VectorDouble a1, double c2, VectorDouble a2, double c3,
			VectorDouble a3, double c4, VectorDouble a4) {
		if ((a1.coordinates.length != a2.coordinates.length) || (a2.coordinates.length != a3.coordinates.length)
				|| (a3.coordinates.length != a4.coordinates.length)) {
			throw new IllegalArgumentException("Input vectors do not have the same length.");
		}
		VectorDouble linearCombo = new VectorDouble(a1.coordinates.length);
		for (int i = 0; i < linearCombo.coordinates.length; i++) {
			linearCombo.coordinates[i] = c1 * a1.coordinates[i] + c2 * a2.coordinates[i] + c3 * a3.coordinates[i] + c4
					* a4.coordinates[i];
		}
		return linearCombo;
	}
	
	/**
	 * Computes a dot product.
	 * 
	 * @param a		vector
	 * return		dot product of the vectors
	 */
	public double dotProduct(VectorDouble a) {
		return dotProduct(this, a);
	}

	/**
	 * Computes a dot product.
	 * 
	 * @param a		first vector
	 * @param b		second vector
	 * return		dot product of the vectors
	 */
	public static double dotProduct(VectorDouble a, VectorDouble b) {
		if (a.coordinates.length != b.coordinates.length) {
			throw new IllegalArgumentException("The input vectors do not have the same length.");
		}

		double result = 0.0;
		for (int i = 0; i < a.coordinates.length; i++) {
			result += a.coordinates[i] * b.coordinates[i];
		}
		return result;
	}

	/**
	 * Computes the square of the norm.
	 * 
	 * return		square of the norm
	 */
	public double normSquared() {
		return this.dotProduct(this);
	}

	/**
	 * Computes the norm.
	 * 
	 * return		norm
	 */
	public double norm() {
		return Math.sqrt(normSquared());
	}

	/**
	 * Makes this vector a unit vector.
	 */
	public void unitVector() {
		double norm = this.norm();
		if (norm == 0) {
			throw new IllegalArgumentException("Norm is zero.");
		}
		this.scalarMult(1 / norm);
	}

	/**
	 * Computes a unit vector.
	 * 
	 * @param a		vector
	 * return		unit vector
	 */
	public static VectorDouble unitVector(VectorDouble a) {
		double norm = a.norm();
		if (norm == 0) {
			throw new IllegalArgumentException("Norm is zero.");
		}
		return scalarMult(1 / norm, a);
	}

	/**
	 * Computes a projection.
	 * 
	 * @param u		defines the projection line
	 */
	public void projection(VectorDouble u) {
		if (u.norm() == 0) {
			throw new IllegalArgumentException("Cannot project onto a vector with norm zero.");
		}
		
		double c = dotProduct(u, this) / u.normSquared();
		for (int i = 0; i < coordinates.length; i++) {
			coordinates[i] = c * u.coordinates[i];
		}
	}
	
	/**
	 * Computes a projection.
	 * 
	 * @param u		defines the projection line
	 * @param v		vector to project
	 * @return		projected vector
	 */
	public static VectorDouble projection(VectorDouble u, VectorDouble v) {
		if (u.norm() == 0) {
			throw new IllegalArgumentException("Cannot project onto a vector with norm zero.");
		}
		return scalarMult(dotProduct(u, v) / u.normSquared(), u);
	}

	/**
	 * Computes a projection.
	 * 
	 * @param u				defines the linear space
	 * @param orthogonal	true if the vectors u are orthogonal
	 */
	public void projection(VectorDouble[] u, boolean orthogonal) {
		VectorDouble[] u_orthogonal = (orthogonal) ? (u) : (doGramSchmidt(u));

		if (u_orthogonal == null) {
			coordinates = new double[coordinates.length];
		} else {
			VectorDouble result = new VectorDouble(coordinates.length);

			for (int i = 0; i < u_orthogonal.length; i++) {
				result.add(VectorDouble.projection(u_orthogonal[i], this));
			}

			coordinates = result.coordinates;
		}
	}

	/**
	 * Computes a projection.
	 * 
	 * @param u		defines the linear space
	 */
	public void projection(VectorDouble[] u) {
		projection(u, false);
	}

	/**
	 * Computes a projection.
	 * 
	 * @param u				defines the linear space
	 * @param v				vector to project
	 * @param orthogonal	true if the vectors u are orthogonal
	 * @return				projected vector
	 */
	public static VectorDouble projection(VectorDouble[] u, VectorDouble v, boolean orthogonal) {
		VectorDouble[] u_orthogonal = (orthogonal) ? (u) : (doGramSchmidt(u));

		if (u_orthogonal == null) {
			return new VectorDouble(v.coordinates.length);
		} else {
			VectorDouble result = new VectorDouble(v.coordinates.length);

			for (int i = 0; i < u_orthogonal.length; i++) {
				result.add(VectorDouble.projection(u_orthogonal[i], v));
			}

			return result;
		}
	}

	/**
	 * Computes a projection.
	 * 
	 * @param u				defines the linear space
	 * @param v				vector to project
	 * @return				projected vector
	 */
	public static VectorDouble projection(VectorDouble[] u, VectorDouble v) {
		return projection(u, v, false);
	}

	/**
	 * Computes a perpendicular projection.
	 * 
	 * @param u		defines the projection line
	 */
	public void perpendicularProjection(VectorDouble u) {
		this.subtract(projection(u, this));
	}

	/**
	 * Computes a perpendicular projection.
	 * 
	 * @param u		defines the projection line
	 * @param v		vector to project
	 * @return		projected vector
	 */
	public static VectorDouble perpendicularProjection(VectorDouble u, VectorDouble v) {
		return subtract(v, projection(u, v));
	}

	/**
	 * Computes a perpendicular projection.
	 * 
	 * @param u				defines the linear space
	 * @param orthogonal	true if the vectors u are orthogonal
	 */
	public void perpendicularProjection(VectorDouble[] u, boolean orthogonal) {
		this.subtract(projection(u, this, orthogonal));
	}

	/**
	 * Computes a perpendicular projection.
	 * 
	 * @param u		defines the linear space
	 */
	public void perpendicularProjection(VectorDouble[] u) {
		perpendicularProjection(u, false);
	}

	/**
	 * Computes a perpendicular projection.
	 * 
	 * @param u				defines the linear space
	 * @param v				vector to project
	 * @param orthogonal	true if the vectors u are orthogonal
	 * @return				projected vector
	 */
	public static VectorDouble perpendicularProjection(VectorDouble[] u, VectorDouble v, boolean orthogonal) {
		return subtract(v, projection(u, v, orthogonal));
	}

	/**
	 * Computes a perpendicular projection.
	 * 
	 * @param u		defines the linear space
	 * @param v		vector to project
	 * @return		projected vector
	 */
	public static VectorDouble perpendicularProjection(VectorDouble u[], VectorDouble v) {
		return perpendicularProjection(u, v, false);
	}
	
	/**
	 * Computes the distance between two vectors.
	 * 
	 * @param a		first vector
	 * @param b		second vector
	 * @return		distance
	 */
	public static double distance(VectorDouble a, VectorDouble b) {
		double tmp, distance = 0;

		if (a.coordinates.length != b.coordinates.length) {
			throw new IllegalArgumentException("The input vectors do not have the same length.");
		}

		for (int i = 0; i < a.coordinates.length; i++) {
			tmp = a.coordinates[i] - b.coordinates[i];
			distance += tmp * tmp;
		}

		return Math.sqrt(distance);
	}
	
	/**
	 * Computes the angle between two vectors.
	 * 
	 * @param a		first vector
	 * @param b		second vector
	 * @return		angle between the two vectors
	 */
	public static double angle(VectorDouble a, VectorDouble b) {
		double normA = a.norm();
		double normB = b.norm();

		if ((normA == 0) || (normB == 0)) {
			throw new IllegalArgumentException("One of the vector norms is zero.");
		}
		
		double cosAngle = VectorDouble.dotProduct(a, b) / (normA * normB);
		if (cosAngle < -1) {
			cosAngle = -1;
		} else if (cosAngle > 1) {
			cosAngle = 1;
		}
		return Math.acos(cosAngle);
	}
	
	/**
	 * Computes the distance from a vector to a line.
	 * 
	 * @param vector	vector
	 * @param a			first point defining line
	 * @param b			second point defining line
	 * @return			distance from the vector to the line
	 */
	public static double distanceToLine(VectorDouble vector, VectorDouble a, VectorDouble b) {
		VectorDouble tmp = VectorDouble.subtract(vector, a);
		tmp.perpendicularProjection((VectorDouble.subtract(b, a)));
		return tmp.norm();
	}
	
	/**
	 * Computes the distance from a vector to a line segment.
	 * 
	 * @param vector	vector
	 * @param a			first point defining line segment
	 * @param b			second point defining line segment
	 * @return			distance from the vector to the line segment
	 */
	public static double distanceToSegment(VectorDouble vector, VectorDouble a, VectorDouble b) {
		VectorDouble vs = subtract(vector, a);
		VectorDouble ve = subtract(vector, b);
		VectorDouble es = subtract(b, a);
		
		if ((dotProduct(vs, es) < 0) && (dotProduct(ve, negate(es))) < 0) {
			return perpendicularProjection(es, vs).norm();
		} else {
			return (vs.norm() < ve.norm()) ? vs.norm() : ve.norm();
		}
	}
	
	/**
	 * Moves the vector to the nearest point on a line segment.
	 * 
	 * @param start			first point defining line segment
	 * @param end			second point defining line segment
	 */
	public void moveToSegment(VectorDouble start, VectorDouble end) {
		VectorDouble a = subtract(this, start);
		VectorDouble b = subtract(this, end);
		VectorDouble c = subtract(end, start);
		
		if ((dotProduct(a, c) < 0) && (dotProduct(b, negate(c))) < 0) {
			a.projection(c);
			a.add(start);
			this.assign(a);
		} else if (a.norm() < b.norm()) {
			this.assign(start);
		} else {
			this.assign(end);
		}
	}
	
	/**
	 * Checks the vectors are of the same dimension and returns the dimension.
	 * 
	 * @param vectors		vectors
	 * @return				consistent dimension
	 */
	public static int checkConsistentDimension(VectorDouble[] vectors) {
		if (vectors == null) {
			throw new IllegalArgumentException("Vectors reference is null.");
		}
		if (vectors.length == 0) {
			throw new IllegalArgumentException("Vectors array is empty.");
		}

		int dimension = 0;
		for (int i = 0; i < vectors.length; i++) {
			if (vectors[i] == null) {
				throw new IllegalArgumentException("Vectors reference is null.");
			}

			if (dimension == 0) {
				dimension = vectors[i].coordinates.length;
			} else {
				if (vectors[i].coordinates.length != dimension) {
					throw new IllegalArgumentException("Dimension is inconsistent.");
				}
			}
		}

		if (dimension == 0) {
			throw new IllegalArgumentException("Dimension is not positive.");
		}
		return dimension;
	}
	
	/**
	 * Performs Gram-Schmidt.
	 * 
	 * @param v		input vectors
	 * @return		vectors which are orthonormal and have the same span as the input vectors
	 */
	public static VectorDouble[] doGramSchmidt(VectorDouble[] v) {
		if (v == null) {
			throw new IllegalArgumentException("The vectors reference is null.");
		}
		int i, j;

		i = 0;
		while (v[i].norm() == 0) {
			if (i < v.length - 1) {
				i++;
			} else {
				return null;
			}
		}

		VectorDouble tmp;
		VectorDouble[] u = new VectorDouble[v.length];
		u[0] = new VectorDouble(v[i]);
		u[0].unitVector();
		int u_size = 1;

		for (i++; i < v.length; i++) {
			tmp = new VectorDouble(v[0].coordinates.length);

			for (j = 0; j < i; j++) {
				tmp.add(projection(u[j], v[i]));
			}

			u[u_size] = subtract(v[i], tmp);

			if ((u[u_size].norm() != 0) && (u_size < u[0].coordinates.length)) {
				u[u_size].unitVector();
				u_size++;
			}
		}

		VectorDouble[] u_new = new VectorDouble[u_size];

		for (i = 0; i < u_size; i++) {
			u_new[i] = u[i];
		}

		return u_new;
	}

	/**
	 * Estimates the tangent space at a node by performing principal component analysis on the incident edges.
	 * 
	 * @param node				node
	 * @param incidentNodes		incident nodes
	 * @param tangentDimension	dimension of the tangent space
	 * @return					orthonormal vectors defining the tangent space at the node
	 */
	public static VectorDouble[] getTangentSpace(VectorDouble node, VectorDouble[] incidentNodes, int tangentDimension) {
		double[][] mArr;
		double[] svalues;
		Matrix m;
		SingularValueDecomposition svd;
		VectorDouble[] vectors;
		VectorDouble tmp;
		int i;

		if (tangentDimension <= 0) {
			throw new IllegalArgumentException("The tangent dimension should be positive.");
		}

		vectors = new VectorDouble[incidentNodes.length];
		for (i = 0; i < incidentNodes.length; i++) {
			vectors[i] = new VectorDouble(incidentNodes[i]);
		}

		if (vectors.length <= tangentDimension) {
			for (i = 0; i < vectors.length; i++) {
				vectors[i].subtract(node);
			}
			return doGramSchmidt(vectors);

		} else {
			tmp = new VectorDouble(node.coordinates.length);
			for (i = 0; i < incidentNodes.length; i++) {
				tmp.add(vectors[i]);
			}

			mArr = new double[vectors.length][];
			tmp.scalarMult(1.0 / vectors.length);

			for (i = 0; i < vectors.length; i++) {
				vectors[i].subtract(tmp);
				mArr[i] = vectors[i].getCoordinates();
			}

			m = new Matrix(mArr, vectors.length, node.coordinates.length);

			svd = m.svd();
			m = svd.getV().transpose();
			mArr = m.getArray();
			svalues = svd.getSingularValues();
			
			for (i = 0; (i < node.coordinates.length) && (i < tangentDimension) && (svalues[i] > 0.0); i++)
				;
			vectors = new VectorDouble[i];

			for (i = 0; i < vectors.length; i++) {
				vectors[i] = new VectorDouble(mArr[i]);
			}

			return vectors;
		}
	}
	
	/**
	 * Converts the vector to a string.
	 * 
	 * @param delimiter		delimiter between coordinates
	 * @param leftBracket	left bracket
	 * @param rightBracket	right bracket
	 * @param format		format for printing coordinates			
	 * @return 				vector as a string
	 */
	public String toString(String delimiter, String leftBracket, String rightBracket, String format) {
		StringWriter s = new StringWriter();
		DecimalFormat dec = new DecimalFormat(format, new DecimalFormatSymbols(Locale.US));
		s.write(leftBracket);

		for (int i = 0; i < this.coordinates.length; i++) {
			s.write(dec.format(coordinates[i]));
			if (i != this.coordinates.length - 1) {
				s.write(delimiter);
			}
		}

		s.write(rightBracket);
		return s.toString();
	}

	/**
	 * Converts the vector to a string, with five decimal places.
	 * 
	 * @param delimiter		delimiter between coordinates
	 * @param leftBracket	left bracket
	 * @param rightBracket	right bracket		
	 * @return 				vector as a string
	 */
	public String toString(String delimiter, String leftBracket, String rightBracket) {
		return toString(delimiter, leftBracket, rightBracket, " 0.00000;-0.00000");
	}

	/**
	 * Converts the vector to a string, with five decimal places and parenthesis brackets.
	 * 
	 * @param delimiter		delimiter between coordinates	
	 * @return 				vector as a string
	 */
	public String toString(String delimiter) {
		return toString(delimiter, "(", ")");
	}
	
	/**
	 * Converts the vector to a string, with five decimal places, parenthesis brackets, and comma delimiter.
	 * 	
	 * @return vector as a string
	 */
	@Override
	public String toString() {
		return toString(", ");
	}

	/**
	 * Converts the vector to a string, with five decimal places.
	 * 
	 * @param comma		true for comma delimiter
	 * @param bracket	true for parenthesis brackets
	 * @return 			vector as a string
	 */
	public String toString(boolean comma, boolean bracket) {
		if (comma) {
			if (bracket) {
				return toString(", ", "(", ")");
			} else {
				return toString(", ", "", "");
			}
		} else {
			if (bracket) {
				return toString(" ", "(", ")");
			} else {
				return toString(" ", "", "");
			}
		}
	}
	
}