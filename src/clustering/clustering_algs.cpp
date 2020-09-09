#include "clustering_algs.h"
#include "../matrix.h"
#include "../random.h"
#include "../union_find.h"
#include <limits>
#include "../IntegralFrechet/IntegralFrechet.h"
#include "../curve_simplification.h"


namespace
{

// TODO: Computes all distances, not only one per pair.
template <typename Comp>
Clustering linkage(Curves const& curves, std::size_t k, int l, Comp comp, distance_t(*dist_func)(Curve, Curve), bool naive_simplification)
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
	while (union_find.getRoots().size() > k) {
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
		auto cl_id = to_cluster_id[union_find.findRoot(curve_id)];
		result[cl_id].curve_ids.push_back(curve_id);
	}

	// We just take the root curves as centers. They don't have any special meaning,
	// but at least we supply some centers.
	for (auto curve_id: union_find.getRoots()) {
		auto cl_id = to_cluster_id[curve_id];
		result[cl_id].center_curve = curves[curve_id].naive_l_simplification(l);
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

Clustering computeClustering(
	Curves const& curves, std::size_t k, int l, ClusterAlg cluster_alg, distance_t(*dist_func)(Curve, Curve), std::string dist_matrix, 
	bool naive_simplification)
{
	switch (cluster_alg) {
	case ClusterAlg::SingleLinkage:
		return singleLinkage(curves, k, l, dist_func, naive_simplification);
	case ClusterAlg::CompleteLinkage:
		return completeLinkage(curves, k, l, dist_func, naive_simplification);
	case ClusterAlg::Gonzalez:
		return runGonzalez(curves, k, l, dist_func, naive_simplification, dist_matrix);
	case ClusterAlg::Pam:
		return pam_with_simplifications(curves, k, l, dist_func, dist_matrix);
	case ClusterAlg::GonzalezPam:
		return gonzalez_pam(curves, k, l, dist_func, dist_matrix);
	}
	
	// ERROR("No matching cluster_alg enum passed.");
}

Clustering singleLinkage(Curves const& curves, std::size_t k, int l, distance_t(*dist_func)(Curve, Curve), bool naive_simplification)
{
	auto min = [](distance_t a, distance_t b) { return std::min<distance_t>(a,b); };
	return linkage(curves, k, l, min, dist_func, naive_simplification);
}

Clustering completeLinkage(Curves const& curves, std::size_t k, int l, distance_t(*dist_func)(Curve, Curve), bool naive_simplification)
{
	auto max = [](distance_t a, distance_t b) { return std::max<distance_t>(a,b); };
	return linkage(curves, k, l, max, dist_func, naive_simplification);
}

Clustering runGonzalez(Curves const& curves, std::size_t k, int l, 
distance_t(*dist_func)(Curve, Curve), bool naive_simplification, std::string dist_matrix_path="")
{
	Clustering result;

	// FrechetLight frechet_light;
	auto max_dist = std::numeric_limits<distance_t>::max();
	std::vector<distance_t> distances_to_center(curves.size(), max_dist);
	ClusterIDs closest_center(curves.size());

	Random random;
	CurveID center_id = random.getUniformInt(0, curves.size()-1);
	std::cout << center_id << std::endl;
	// CurveID center_id = random.getUniformInt<Curves::size_type>(0, curves.size() - 1);

	// load the distance matrix if provided
	CurveSimpMatrix dist_matrix = dist_matrix_path == "" ? CurveSimpMatrix()
	: CurveSimpMatrix(dist_matrix_path);


	// add as center and update closest distances to center
	// auto center_curve = simplify(curves[center_id], l);
	auto center_curve = naive_simplification ? curves[center_id].naive_l_simplification(l) : simplify(curves[center_id], l, dist_func);
	
	result.push_back({{}, center_curve, std::numeric_limits<distance_t>::max(), center_id});
	for (CurveID curve_id = 0; curve_id < curves.size(); ++curve_id) {
		
		auto& current_dist = distances_to_center[curve_id];
		// auto new_dist = compute_integral_frechet_distance(center_curve, curves[curve_id]);
		auto new_dist = dist_matrix_path == "" ? dist_func(center_curve, curves[curve_id]) 
		: dist_matrix.at(curve_id, center_id);
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
	
	while (result.size() < k) {
		auto center_it = std::max_element(distances_to_center.begin(), distances_to_center.end());
		auto cid = static_cast<std::size_t>(std::distance(distances_to_center.begin(), center_it));
		
		Curve cent_curve = naive_simplification ? 
		curves[cid].naive_l_simplification(l) : 
		simplify(curves[cid], l, dist_func);
		

		result.push_back({{}, cent_curve, std::numeric_limits<distance_t>::max(), cid});
		distances_to_center[cid] = 0;
		closest_center[cid] = result.size() - 1;
		
		for (CurveID curve_id = 0; curve_id < curves.size(); ++curve_id) {

			auto& current_dist = distances_to_center[curve_id];

			auto new_dist = dist_matrix_path == "" ? dist_func(cent_curve, curves[curve_id]) 
			: dist_matrix.at(curve_id, cid);
			if (new_dist < current_dist) {
				current_dist = new_dist;
				closest_center[curve_id] = result.size()-1;
			}
		}
	}

	for (CurveID curve_id = 0; curve_id < curves.size(); ++curve_id) {
		auto cluster_id = closest_center[curve_id];
		result[cluster_id].curve_ids.push_back(curve_id);
	}

	return result;
}

void updateClustering(Curves const& curves, Clustering& clustering, distance_t(*dist_func)(Curve, Curve), CurveSimpMatrix* distance_matrix)
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
			auto new_dist = distance_matrix == nullptr ? dist_func(curves[curve_id], center_curve) :
			distance_matrix->at(curve_id, clustering[cluster_id].center_id);

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

distance_t calcDiameter(Curves const& curves, CurveIDs const& /* curve_ids */, distance_t(*dist_func)(Curve, Curve))
{
	// FrechetLight frechet_light;
	distance_t max_distance = 0.;
	for (CurveID curve_id1 = 0; curve_id1 < curves.size(); ++curve_id1) {
		for (CurveID curve_id2 = curve_id1 + 1; curve_id2 < curves.size(); ++curve_id2) {
			// auto dist = frechet_light.calcDistance(curves[curve_id1], curves[curve_id2]);
			// auto dist = compute_integral_frechet_distance(curves[curve_id1], curves[curve_id2]);
			auto dist = dist_func(curves[curve_id1], curves[curve_id2]);
			max_distance = std::max(max_distance, dist);
		}
	}

	return max_distance;
}

Clustering pam_with_simplifications(Curves const& curves, std::size_t k, int l, distance_t(*dist_func)(Curve, Curve), std::string matrix_file_name) {
	Curves simplifications = Curves();
	for (auto curve: curves) {
		simplifications.push_back(simplify(curve, l, dist_func));
	}

	CurveSimpMatrix distance_matrix = matrix_file_name == "" ? 
	CurveSimpMatrix(curves, simplifications, dist_func) :
	CurveSimpMatrix(matrix_file_name);

	std::vector<size_t> medoids;

	std::vector<size_t> initial_centers = clustering::pam_simp::compute(curves.size(), k, distance_matrix, medoids);

	Clustering clustering = Clustering();

	for (auto center_id: initial_centers) {
		// std::cout << curves[center_id].name() << "\n";
		clustering.push_back({{}, simplify(curves[center_id], l, dist_func), std::numeric_limits<distance_t>::max(), center_id});
	}

	updateClustering(curves, clustering, dist_func, matrix_file_name == "" ? nullptr : &distance_matrix);

	return clustering;
}

Clustering pam_with_centering(Curves const& curves, std::size_t k, int l, distance_t(*dist_func)(Curve, Curve), std::string matrix_file_name) {
	
	Curves simplifications = Curves();
	for (auto curve: curves) {
		simplifications.push_back(curve.naive_l_simplification(l));
	}

	// std::cout << "computed simplificactions... \n";

	CurveSimpMatrix distance_matrix = matrix_file_name == "" ? 
	CurveSimpMatrix(curves, simplifications, dist_func) :
	CurveSimpMatrix(matrix_file_name);

	// std::cout << "read distance matrix \n";

	std::vector<size_t> medoids;

	std::vector<size_t> initial_centers = clustering::pam_simp::compute(curves.size(), k, distance_matrix, medoids);

	// std::cout << "initial centers: " << initial_centers.size() << "\n";
	for (auto id: initial_centers) {
		std::cout << id << "\n";
	}

	Clustering clustering = Clustering();

	for (auto center_id: initial_centers) {
		std::cout << curves[center_id].name() << "\n";
		clustering.push_back({{}, simplifications[center_id], std::numeric_limits<distance_t>::max()});
	}

	updateClustering(curves, clustering, dist_func);

	for (auto cluster: clustering) {
		// std::cout << "cluster size: " << cluster.curve_ids.size() << "\n";
	}

	distance_t cost_after_pam = 0;
	for (auto& cluster: clustering) {
		cluster.cost = calcC2CDist(curves, cluster.center_curve, cluster.curve_ids, C2CDist::Median, dist_func);
		cost_after_pam += cluster.cost;
	}

	std::cout << "cost after pam: " << cost_after_pam << "\n";

	std::cout << clustering[0].cost << "\n";

	calcFSACenters(curves, clustering, l, dist_func, C2CDist::Median, CenterCurveUpdateMethod::frechetMean);

	std::cout << "updated centers\n";

	updateClustering(curves, clustering, dist_func);


	distance_t cost_after_centering = 0;
	for (auto& cluster: clustering) {
		cost_after_centering += calcC2CDist(curves, cluster.center_curve, cluster.curve_ids, C2CDist::Median, dist_func);
	}

	std::cout << "cost after centering: " << cost_after_centering << "\n";

	return clustering;
}

distance_t kMedianCost(Curves const& curves, Clustering const& clustering, distance_t(*dist_func)(Curve, Curve), CurveSimpMatrix& dist_matrix) {
	distance_t cost = 0;

	for (auto& cluster: clustering) {
		cost += calcC2CDist(curves, cluster.center_curve, cluster.curve_ids, C2CDist::Median, dist_func);
	}

	return cost;
}

distance_t kMedianCostMat(Curves const& curves, Clustering const& clustering, CurveSimpMatrix& dist_matrix) {

	distance_t cost = 0;

	for (auto cluster: clustering) {
		for (auto curve_id: cluster.curve_ids) {
			cost += dist_matrix.at(curve_id, cluster.center_id);
		}
	}

	return cost;
}

Clustering gonzalez_pam(Curves const& curves, std::size_t k, int l, distance_t(*dist_func)(Curve, Curve), std::string matrix_file_name) {
	Clustering initial_clustering = runGonzalez(curves, k, l, dist_func, false, matrix_file_name);
	// load the distance matrix if provided
	std::vector<size_t> initial_centers = std::vector<size_t>();

	// load the distance matrix
	CurveSimpMatrix dist_matrix = matrix_file_name == "" ? CurveSimpMatrix()
	: CurveSimpMatrix(matrix_file_name);

	for (auto cluster: initial_clustering) {
		initial_centers.push_back(cluster.center_id);
	}

	std::vector<size_t> new_centers = clustering::pam_simp::compute(
		curves.size(), k, dist_matrix, initial_centers
	);

	Clustering clustering = Clustering();

	for (auto center_id: new_centers) {
		clustering.push_back({
			{}, simplify(curves[center_id], l, dist_func), std::numeric_limits<distance_t>::max(), center_id
		});
	}

	updateClustering(curves, clustering, dist_func, matrix_file_name == "" ? nullptr : &dist_matrix);
	
	return clustering;
}
