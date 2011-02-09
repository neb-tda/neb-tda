package edu.stanford.math.nebtda;

import java.util.*;

/**
 * Abstract class for a clustering method.
 * 
 * @author Henry Adams
 * @author Atanas Atanasov
 */
public abstract class Clustering<T> {
	
	protected MetricSpace<T> metric;

	/**
	 * Constructor.
	 * 
	 * @param metric		metric
	 */
	public Clustering(MetricSpace<T> metric) {
		setMetric(metric);
	}

	/**
	 * Sets the metric.
	 * 
	 * @param metric		metric
	 */
	public void setMetric(MetricSpace<T> metric) {
		if (metric == null) {
			throw new IllegalArgumentException("The metric reference is null.");
		}
		this.metric = metric;
	}

	/**
	 * Gets the metric.
	 * 
	 * @return metric
	 */
	public MetricSpace<T> getMetric() {
		return metric;
	}
	
	/**
	 * Clusters a set of points.
	 * 
	 * @param points	points to cluster
	 * @return 			a collection of clusters
	 */
	public abstract Collection<Cluster<T>> cluster(T[] points);

}
