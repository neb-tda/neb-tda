package edu.stanford.math.nebtda;

/**
 * Abstract class for choosing a representative from a cluster.
 * 
 * @author Henry Adams
 * @author Atanas Atanasov
 */
public abstract class ClustRep<T> {

	/**
	 * Empty constructor.
	 */
	public ClustRep() {
	}

	/**
	 * Gets the representative of the cluster.
	 * 
	 * @param cluster		cluster
	 * @return				representative
	 */
	public abstract T getRepresentative(Cluster<T> cluster);

}
