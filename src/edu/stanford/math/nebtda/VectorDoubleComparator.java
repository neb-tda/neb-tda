package edu.stanford.math.nebtda;

import java.util.Comparator;

/**
 * Compares vectors, based on the norm of the coordinate at a specified index.
 * Currently used only by the method OpticalImagePatches() in class Main.
 * 
 * @author Henry Adams
 * @author Atanas Atanasov
 */
public class VectorDoubleComparator implements Comparator<VectorDouble> {
	
	protected int index;
	
	/**
	 * Constructor.
	 * 
	 * @param index		the index used to compare vectors
	 */
	public VectorDoubleComparator(int index) {
		this.index = index;
	}
	
	/**
	 * Returns 1 if the norm of the specified index of vector a is less than that of vector b.
	 * Returns 0 if the norms are equal. Else returns -1.
	 * 
	 * @param a		first vector
	 * @param b		second vector
	 * @return		1, 0, or -1
	 */
	public int compare(VectorDouble a, VectorDouble b) {
		double value_a = Math.abs(a.coordinates[index]);
		double value_b = Math.abs(b.coordinates[index]);
		if (value_a < value_b) {
			return 1;
		} else if (value_a > value_b) {
			return -1;
		} else {
			return 0;
		}
		
	}

}
