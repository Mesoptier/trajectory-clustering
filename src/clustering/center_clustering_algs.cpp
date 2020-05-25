#include "center_clustering_algs.h"


namespace
{

Clustering computeCenterClusteringRound(
	Curves const& curves, int k, int l, ClusterAlg cluster_alg, CenterAlg center_alg,
	distance_t(*initial_dist)(Curve, Curve), distance_t(*dist_func)(Curve, Curve), std::string dist_matrix)
{
	auto clustering = computeClustering(curves, k, l, ClusterAlg::Gonzalez, initial_dist, dist_matrix, false);
	updateClustering(curves, clustering, dist_func);

	// iterate as long as there are new centers
	int count = 1;
	int const max_count = 1;
	while (count <= max_count && computerCenters(curves, clustering, l, center_alg, dist_func)) {
		updateClustering(curves, clustering, dist_func);
		++count;
	}
	// std::cout << "Number of iterations: " << count << std::endl;

	return clustering;
}

} // end anonymous namespace

Clustering computeCenterClustering(
	Curves const& curves, int k, int l, ClusterAlg cluster_alg, CenterAlg center_alg, distance_t(*initial_dist)(Curve, Curve), distance_t(*dist_func)(Curve, Curve), std::string dist_matrix, int max_rounds)
{
	Clustering min_clustering;
	distance_t min_cost = std::numeric_limits<distance_t>::max();

	for (int round = 0; round < max_rounds; ++round) {
		auto clustering = computeCenterClusteringRound(curves, k, l, cluster_alg, center_alg, initial_dist, dist_func, dist_matrix);
		distance_t cost_sum = 0.;
		for (auto const& cluster: clustering) {
			cost_sum += cluster.cost;
		}
		if (cost_sum < min_cost) {
			min_clustering = std::move(clustering);
			min_cost = cost_sum;
		}
	}

	// remove empty clusters
	for (ClusterID cluster_id = 0; cluster_id < min_clustering.size(); ++cluster_id) {
		if (min_clustering[cluster_id].curve_ids.empty()) {
			std::swap(min_clustering[cluster_id], min_clustering.back());
			min_clustering.pop_back();
		}
	}
	
	return min_clustering;
}
