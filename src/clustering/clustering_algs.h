#pragma once

#include "../basic_types.h"
#include "../Curve.h"
#include "center_algs.h"
#include "../CurveSimpMatrix.h"
#include "pam_with_simplifications.h"

// #include "../curve_simplification.h"
// #include "../compute_distance.h"

using Curves = std::vector<Curve>;

enum class ClusterAlg {
	SingleLinkage,
	CompleteLinkage,
	Gonzalez,
	Pam,
	GonzalezPam
};
std::string toString(ClusterAlg cluster_alg);

Clustering computeClustering(
	Curves const& curves, std::size_t k, int l, ClusterAlg cluster_alg, distance_t(*dist_func)(Curve, Curve), std::string dist_matrix, 
	bool naive_simplification=false);

Clustering singleLinkage(Curves const& curves, std::size_t k, int l, distance_t(*dist_func)(Curve, Curve), bool naive_simplification);
Clustering completeLinkage(Curves const& curves, std::size_t k, int l, distance_t(*dist_func)(Curve, Curve), bool naive_simplification);
Clustering runGonzalez(Curves const& curves, std::size_t k, int l, 
distance_t(*dist_func)(Curve, Curve), bool naive_simplification, std::string dist_matrix_path);

// assign curves to closest clusters
void updateClustering(Curves const& curves, Clustering& clustering, distance_t(*dist_func)(Curve, Curve), CurveSimpMatrix* distance_matrix=NULL);

distance_t calcDiameter(Curves const& curves, CurveIDs const& curve_ids, distance_t(*dist_func)(Curve, Curve));

distance_t kMedianCost(Curves const& curves, Clustering const& clustering, distance_t(*dist_funct)(Curve, Curve));

distance_t kMedianCostMat(Curves const& curves, Clustering const& clustering, CurveSimpMatrix& dist_matrix);

Clustering pam_with_simplifications(Curves const& curves, std::size_t k, int l, distance_t(*dist_func)(Curve, Curve), std::string matrix_file_name="");

Clustering pam_with_centering(Curves const& curves, std::size_t k, int l, distance_t(*dist_func)(Curve, Curve), std::string matrix_file_name);

Clustering gonzalez_pam(Curves const& curves, std::size_t k, int l, distance_t(*dist_func)(Curve, Curve), std::string matrix_file_name="");