#include "clustering_algs.h"
#include "../matrix.h"
#include "../random.h"
#include "../union_find.h"
#include <limits>
#include "../IntegralFrechet/IntegralFrechet.h"
// #include "../curve_simplification.h"

namespace
{


distance_t compute_integral_frechet_distance(Curve curve_1, Curve curve_2) {
	/*const auto result_alt = IntegralFrechet(curve_1.simplify(true), curve_2, ParamMetric::LInfinity_NoShortcuts, 10, nullptr).compute_matching();
    const auto band = MatchingBand(curve_1, curve_2, result_alt.matching, 1);
	const auto result = IntegralFrechet(curve_1, curve_1, ParamMetric::LInfinity_NoShortcuts, 1, &band).compute_matching();
	return result.cost;*/
	
	distance_t cost = IntegralFrechet(
        curve_1, curve_2, ParamMetric::LInfinity_NoShortcuts, 1, nullptr
    ).compute_matching()
    .cost;
    return cost;
}


// TODO: Computes all distances, not only one per pair.
template <typename Comp>
Clustering linkage(Curves const& curves, int k, int l, Comp comp, distance_t(*dist_func)(Curve, Curve))
{
	// compute all pairwise Fréchet distances
	// FrechetLight frechet_light;
	Matrix<distance_t> dist_matrix(curves.size(), curves.size());
	for (CurveID curve_id1 = 0; curve_id1 < curves.size(); ++curve_id1) {
		for (CurveID curve_id2 = 0; curve_id2 < curves.size(); ++curve_id2) {
			if (curve_id1 == curve_id2) {
				dist_matrix(curve_id1, curve_id2) = 0;
			}
			else if (curve_id1 < curve_id2) {
				// auto dist = frechet_light.calcDistance(curves[curve_id1], curves[curve_id2]);
				// auto dist = compute_integral_frechet_distance(curves[curve_id1], curves[curve_id2]);
				auto dist = dist_func(curves[curve_id1], curves[curve_id2]);
				dist_matrix(curve_id1,curve_id2) = dist;
				dist_matrix(curve_id2,curve_id1) = dist;
			}
		}
	}

	// create initial clusters
	CurveIDs base_set(curves.size());
	std::iota(base_set.begin(), base_set.end(), 0);
	UnionFind<CurveID> union_find(base_set);

	// merge clusters until there are exactly k
	while ((int)union_find.getRoots().size() > k) {
		// find two clusters to merge
		distance_t min_dist = std::numeric_limits<distance_t>::max();
		CurveID min_id1 = CurveID(), min_id2 = CurveID();
		for (auto curve_id1: union_find.getRoots()) {
			for (auto curve_id2: union_find.getRoots()) {
				if (curve_id1 == curve_id2) { continue; }

				auto new_dist = dist_matrix(curve_id1, curve_id2);
				if (new_dist < min_dist) {
					min_dist = new_dist;
					min_id1 = curve_id1;
					min_id2 = curve_id2;
				}
			}
		}

		// merge clusters and adapt distances
		min_id1 = union_find.findRoot(min_id1);
		min_id2 = union_find.findRoot(min_id2);
		auto new_root_id = union_find.uniteSets(min_id1, min_id2);
		auto new_child_id = (new_root_id == min_id1 ? min_id2 : min_id1);

		for (auto curve_id: union_find.getRoots()) {
			auto root_dist = dist_matrix(new_root_id, curve_id);
			auto child_dist = dist_matrix(new_child_id, curve_id);
			dist_matrix(new_root_id, curve_id) = comp(root_dist, child_dist);
		}
	}

	// construct the result
	Clustering result(union_find.getRoots().size());

	std::unordered_map<CurveID, std::size_t> to_cluster_id;
	ClusterID cluster_id = 0;
	for (auto curve_id: union_find.getRoots()) {
		to_cluster_id[curve_id] = cluster_id;
		++cluster_id;
	}
	for (CurveID curve_id = 0; curve_id < curves.size(); ++curve_id) {
		auto cluster_id = to_cluster_id[union_find.findRoot(curve_id)];
		result[cluster_id].curve_ids.push_back(curve_id);
	}

	// We just take the root curves as centers. They don't have any special meaning,
	// but at least we supply some centers.
	for (auto curve_id: union_find.getRoots()) {
		auto cluster_id = to_cluster_id[curve_id];
		result[cluster_id].center_curve = curves[curve_id].simplify(true);
		// result[cluster_id].center_curve = simplify(curves[curve_id], l);
	}

	return result;
}

} // end anonymous namespace

std::string toString(ClusterAlg cluster_alg) {
	switch(cluster_alg) {
	case ClusterAlg::SingleLinkage: return "SingleLinkage";
	case ClusterAlg::CompleteLinkage: return "CompleteLinkage";
	case ClusterAlg::Gonzalez: return "Gonzalez";
	}

	// ERROR("Unknown cluster_alg.");
}

Clustering computeClustering(Curves const& curves, int k, int l, ClusterAlg cluster_alg, distance_t(*dist_func)(Curve, Curve))
{
	switch (cluster_alg) {
	case ClusterAlg::SingleLinkage:
		return singleLinkage(curves, k, l, dist_func);
	case ClusterAlg::CompleteLinkage:
		return completeLinkage(curves, k, l, dist_func);
	case ClusterAlg::Gonzalez:
		return runGonzalez(curves, k, l, dist_func);
	}

	
	// ERROR("No matching cluster_alg enum passed.");
}

Clustering singleLinkage(Curves const& curves, int k, int l, distance_t(*dist_func)(Curve, Curve))
{
	auto min = [](distance_t a, distance_t b) { return std::min<distance_t>(a,b); };
	return linkage(curves, k, l, min, dist_func);
}

Clustering completeLinkage(Curves const& curves, int k, int l, distance_t(*dist_func)(Curve, Curve))
{
	auto max = [](distance_t a, distance_t b) { return std::max<distance_t>(a,b); };
	return linkage(curves, k, l, max, dist_func);
}

Clustering runGonzalez(Curves const& curves, int k, int l, distance_t(*dist_func)(Curve, Curve))
{
	Clustering result;

	// FrechetLight frechet_light;
	auto max_dist = std::numeric_limits<distance_t>::max();
	std::vector<distance_t> distances_to_center(curves.size(), max_dist);
	ClusterIDs closest_center(curves.size());

	Random random;
	CurveID center_id = random.getUniformInt(0, curves.size()-1);

	// add as center and update closest distances to center
	// auto center_curve = simplify(curves[center_id], l);
	auto center_curve = curves[center_id].simplify(true);
	
	result.push_back({{}, center_curve});
	for (CurveID curve_id = 0; curve_id < curves.size(); ++curve_id) {
		auto& current_dist = distances_to_center[curve_id];
		// auto new_dist = compute_integral_frechet_distance(center_curve, curves[curve_id]);
		auto new_dist = dist_func(center_curve, curves[curve_id]);
		if (new_dist < current_dist) {
			current_dist = new_dist;
			closest_center[curve_id] = result.size()-1;
		}

		/*
		if (frechet_light.lessThanWithFilters(current_dist, center_curve, curves[curve_id])) {
			current_dist = frechet_light.calcDistance(center_curve, curves[curve_id]);
			closest_center[curve_id] = result.size()-1;
		}*/
	}

	while ((int)result.size() < k) {
		auto center_it = std::max_element(distances_to_center.begin(), distances_to_center.end());
		auto center_id = std::distance(distances_to_center.begin(), center_it);
		// auto center_curve = simplify(curves[center_id], l);
		Curve center_curve = curves[center_id].simplify(true);
		result.push_back({{}, center_curve});
		for (CurveID curve_id = 0; curve_id < curves.size(); ++curve_id) {
			auto& current_dist = distances_to_center[curve_id];

			// auto new_dist = compute_integral_frechet_distance(center_curve, curves[curve_id]);
			auto new_dist = dist_func(center_curve, curves[curve_id]);
			if (new_dist < current_dist) {
				current_dist = new_dist;
				closest_center[curve_id] = result.size()-1;
			}

			/*
			if (frechet_light.lessThanWithFilters(current_dist, center_curve, curves[curve_id])) {
				current_dist = frechet_light.calcDistance(center_curve, curves[curve_id]);
				closest_center[curve_id] = result.size()-1;
			}*/
		}
	}

	for (CurveID curve_id = 0; curve_id < curves.size(); ++curve_id) {
		auto cluster_id = closest_center[curve_id];
		result[cluster_id].curve_ids.push_back(curve_id);
	}

	return result;
}

void updateClustering(Curves const& curves, Clustering& clustering, distance_t(*dist_func)(Curve, Curve))
{
	// FrechetLight frechet_light;

	// clear clusters
	for (auto& cluster: clustering) {
		cluster.curve_ids.clear();
	}

	// compute new clusters
	for (CurveID curve_id = 0; curve_id < curves.size(); ++curve_id) {
		distance_t min_dist = std::numeric_limits<distance_t>::max();
		ClusterID min_cluster_id = ClusterID();
		for (ClusterID cluster_id = 0; cluster_id < clustering.size(); ++cluster_id) {
			auto const& center_curve = clustering[cluster_id].center_curve;
			
			// auto new_dist = compute_integral_frechet_distance(curves[curve_id], center_curve);
			auto new_dist = dist_func(curves[curve_id], center_curve);
			if (new_dist < min_dist) {
				min_dist = new_dist;
				min_cluster_id = cluster_id;
			}
			
			/*if (frechet_light.lessThanWithFilters(min_dist, curves[curve_id], center_curve)) {
				min_dist = frechet_light.calcDistance(curves[curve_id], center_curve);
				min_cluster_id = cluster_id;
			}*/
		}

		clustering[min_cluster_id].curve_ids.push_back(curve_id);
	}
}

distance_t calcDiameter(Curves const& curves, CurveIDs const& curve_ids, distance_t(*dist_func)(Curve, Curve))
{
	// FrechetLight frechet_light;
	distance_t max_distance = 0.;
	for (CurveID curve_id1 = 0; curve_id1 < curves.size(); ++curve_id1) {
		for (CurveID curve_id2 = 0; curve_id2 < curves.size(); ++curve_id2) {
			if (curve_id1 == curve_id2) {
				continue;
			}
			else if (curve_id1 < curve_id2) {
				// auto dist = frechet_light.calcDistance(curves[curve_id1], curves[curve_id2]);
				// auto dist = compute_integral_frechet_distance(curves[curve_id1], curves[curve_id2]);
				auto dist = dist_func(curves[curve_id1], curves[curve_id2]);
				max_distance = std::max(max_distance, dist);
			}
		}
	}

	return max_distance;
}

Clustering pam_with_centering(Curves const& curves, int k, int l, distance_t(*dist_func)(Curve, Curve)) {
	
	Curves simplifications = Curves();
	for (auto curve: curves) {
		simplifications.push_back(curve.simplify(true));
	}

	std::cout << "computed simplificactions... \n";

	CurveSimpMatrix distance_matrix = CurveSimpMatrix(curves, simplifications, false);

	std::vector<size_t> initial_centers = clustering::pam_simp::compute(curves.size(), k, distance_matrix);

	std::cout << "initial centers: " << initial_centers.size() << "\n";
	for (auto id: initial_centers) {
		std::cout << id << "\n";
	}

	Clustering clustering = Clustering();

	for (auto center_id: initial_centers) {
		std::cout << curves[center_id].name() << "\n";
		clustering.push_back({{}, simplifications[center_id], 0});
	}

	updateClustering(curves, clustering, dist_func);

	for (auto cluster: clustering) {
		std::cout << "cluster size: " << cluster.curve_ids.size() << "\n";
	}

	distance_t cost_after_pam = 0;
	for (auto& cluster: clustering) {
		cluster.cost = calcC2CDist(curves, cluster.center_curve, cluster.curve_ids, C2CDist::Median);
		cost_after_pam += cluster.cost;
	}

	std::cout << "cost after pam: " << cost_after_pam << "\n";

	std::cout << clustering[0].cost << "\n";

	calcFSACenters(curves, clustering, l, C2CDist::Median);

	std::cout << "updated centers\n";

	for (auto cluster: clustering) {
		updateClustering(curves, clustering, dist_func);
	}

	distance_t cost_after_centering = 0;
	for (auto& cluster: clustering) {
		cost_after_centering += calcC2CDist(curves, cluster.center_curve, cluster.curve_ids, C2CDist::Median);
	}

	std::cout << "cost after centering: " << cost_after_centering << "\n";

	return clustering;
}