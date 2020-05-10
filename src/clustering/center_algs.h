#pragma once

#include "../basic_types.h"
#include "../defs.h"
#include "../DTW/dtw.h"

#include <string>

using Curves = std::vector<Curve>;

enum class CenterAlg {
	kMedian,
	kMeans,
	kCenter,
	fCenter,
	fMean,
	dtwMean,
	avFCenter,
	newCenterUpdate,
	newCenterUpdate2,
	naiveCenterUpdate,
	ensembleMethod1
};

enum class CenterCurveUpdateMethod {
	frechetCentering,
	frechetMean,
	dtwMean,
	avFCenter,
	newCenterUpdate,
	newCenterUpdate2,
	naiveCenterUpdate,
	ensembleMethod1
};

std::string toString(CenterAlg center_alg);

enum class C2CDist {
	Median,
	Mean,
	Max,
};

distance_t calcC2CDist(Curves const& curves, Curve const& center_curve, CurveIDs const& curve_ids, C2CDist c2c_dist, distance_t(*dist_func)(Curve, Curve));

bool computerCenters(Curves const& curves, Clustering& clustering, int l, CenterAlg center_alg, distance_t(*dist_func)(Curve, Curve));

bool calcKMedianCenters(Curves const& curves, Clustering& clustering, int l, distance_t(*dist_func)(Curve, Curve));
bool calcKMeansCenters(Curves const& curves, Clustering& clustering, int l, distance_t(*dist_func)(Curve, Curve));
bool calcKCenterCenters(Curves const& curves, Clustering& clustering, int l, distance_t(*dist_func)(Curve, Curve));
bool calcFSACenters(
	Curves const& curves, Clustering& clustering, int l, distance_t(*dist_func)(Curve, Curve), C2CDist cluster_dist = C2CDist::Max, CenterCurveUpdateMethod method=CenterCurveUpdateMethod::frechetMean);
bool naiveCenterUpdate(Curves const& curves, Clustering& clustering, distance_t(*dist_func)(Curve, Curve), C2CDist c2c_dist);

// First computes a new center via the naiveCenterUpdateMethod
// the frechetMean update method is then applied
bool ensembleMethod1(Curves const& curves, Clustering& clustering, distance_t(*dist_func)(Curve, Curve), C2CDist c2c_dist);
Points matching_of_vertices(Curve curve_1, Curve curve_2);
