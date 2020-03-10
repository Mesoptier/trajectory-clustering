#pragma once

#include "../basic_types.h"
#include "../defs.h"

#include <string>

using Curves = std::vector<Curve>;

enum class CenterAlg {
	kMedian,
	kMeans,
	kCenter,
	FSA,
};
std::string toString(CenterAlg center_alg);

enum class C2CDist {
	Median,
	Mean,
	Max,
};

distance_t calcC2CDist(Curves const& curves, Curve const& center_curve, CurveIDs const& curve_ids, C2CDist c2c_dist);

// bool computerCenters(Curves const& curves, Clustering& clustering, int l, CenterAlg center_alg);

// bool calcKMedianCenters(Curves const& curves, Clustering& clustering, int l);
// bool calcKMeansCenters(Curves const& curves, Clustering& clustering, int l);
// bool calcKCenterCenters(Curves const& curves, Clustering& clustering, int l);
bool calcFSACenters(
	Curves const& curves, Clustering& clustering, int l, C2CDist cluster_dist = C2CDist::Max);

Points matching_of_vertices(Curve curve_1, Curve curve_2);
