package edu.stanford.math.nebtda;

/**
 * Abstract class for a cell.
 * 
 * @author Henry Adams
 * @author Atanas Atanasov
 */
public abstract class Cell implements Value<Cell> {
	
	protected double value;
	
	/**
	 * Sets the density value.
	 * 
	 * @param value 	density value
	 */
	public void setValue(double value) {
		this.value = value;
	}
	
	/**
	 * Gets the density value.
	 * 
	 * @return density value
	 */
	public double getValue() {
		return value;
	}
	
	/**
	 * Gets the density value.
	 * 
	 * @return density value
	 */
	public double getValue(Cell cell) {
		return cell.value;
	}
	
	/**
	 * Computes the density value.
	 * 
	 * @param density	density estimator
	 * @return 			density value
	 */
	public abstract double computeValue(Density density);

}
