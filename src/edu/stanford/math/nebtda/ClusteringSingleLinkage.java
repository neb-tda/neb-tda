package edu.stanford.math.nebtda;

import java.util.*;

/**
 * Implements single-linkage clustering.
 * 
 * @author Henry Adams
 * @author Atanas Atanasov
 */
public class ClusteringSingleLinkage<T> extends Clustering<T> {

	protected double clusteringThreshold;

	/**
	 * Constructor.
	 * 
	 * @param metric		metric
	 * @param epsilon		epsilon
	 */
	public ClusteringSingleLinkage(MetricSpace<T> metric, double epsilon) {
		super(metric);
		setClusteringThreshold(epsilon);
	}

	/**
	 * Sets the clustering threshold.
	 * 
	 * @param clusteringThreshold		clustering threshold
	 */
	public void setClusteringThreshold(double clusteringThreshold) {
		if (clusteringThreshold <= 0) {
			throw new IllegalArgumentException("Epsilon is not positive.");
		}
		this.clusteringThreshold = clusteringThreshold;
	}

	/**
	 * Gets the clustering threshold.
	 * 
	 * @return clustering threshold
	 */
	public double getClusteringThreshold() {
		return clusteringThreshold;
	}

	/**
	 * Performs single-linkage clustering.
	 * 
	 * @param points 	the points to cluster
	 * @return 			collection of clusters
	 */
	@Override
	public Collection<Cluster<T>> cluster(T[] points) {
		LinkedList<Cluster<T>> result = new LinkedList<Cluster<T>>();
		ListIterator<Cluster<T>> it;
		Cluster<T> firstCluster, cluster;
		int i;

		for (i = 0; i < points.length; i++) {
			firstCluster = null;

			for (it = result.listIterator(); it.hasNext();) {
				cluster = it.next();

				for (T t : cluster) {
					if (metric.distance(points[i], t) < clusteringThreshold) {
						if (firstCluster == null) {
							firstCluster = cluster;
							cluster.add(points[i]);
						} else {
							firstCluster.addAll(cluster);
							it.remove();
						}
						break;
					}
				}
			}

			if (firstCluster == null) {
				result.add(new Cluster<T>(points[i]));
			}
		}

		return result;
	}

}
