package edu.stanford.math.nebtda;

import java.util.LinkedList;

/**
 * Implements a cluster.
 * 
 * @author Henry Adams
 * @author Atanas Atanasov
 */
public class Cluster<T> extends LinkedList<T> {
	
	private static final long serialVersionUID = 2210044184293026628L;

	/**
	 * Constructor.
	 */
	public Cluster() {
		super();
	}

	/**
	 * Constructor which adds a single point.
	 * 
	 * @param point		single point in new cluster
	 */
	public Cluster(T point) {
		this();
		this.add(point);
	}

//	/**
//	 * Constructor which copies a cluster.
//	 * 
//	 * @param cluster		cluster to copy
//	 */
//	public Cluster(Cluster<T> cluster) {
//		super((Collection<T>) cluster);
//	}

//	/**
//	 * Constructor which combines multiple clusters.
//	 * 
//	 * @param collections		collection of clusters to combine
//	 */
//	public Cluster(Collection<Cluster<T>> collections) {
//		this();
//		for (Collection<T> c : collections) {
//			this.addAll(c);
//		}
//	}

}
