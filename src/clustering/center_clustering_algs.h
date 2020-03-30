#pragma once

#include "center_algs.h"
#include "clustering_algs.h"
#include "../Curve.h"

Clustering computeCenterClustering(
	Curves const& curves, int k, int l, ClusterAlg cluster_alg, CenterAlg center_alg, distance_t(*dist_func)(Curve, Curve), std::string dist_matrix, 
	int max_rounds = 10);
