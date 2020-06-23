#pragma once

#include "center_update.h"
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
	ensembleMethod1,
	dba,
	cdba,
	wedge,
	wedge_2,
	regression,
	regression_3d,
	cdbaChar,
	cdbaPigeon,
	dbaChar,
	dbaPigeon,
	wedgeChar,
	wedgePigeon
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
bool dba(Curves const& curves, Clustering& clustering, distance_t(*dist_func)(Curve, Curve), C2CDist c2c_dist);

template <int res>
bool cdba(Curves const& curves, Clustering& clustering, distance_t(*dist_func)(Curve, Curve), C2CDist c2c_dist);

bool ensembleMethod1(Curves const& curves, Clustering& clustering, distance_t(*dist_func)(Curve, Curve), C2CDist c2c_dist);
bool updateCenters(Curves const& curves, Clustering& clustering, distance_t(*dist_func)(Curve, Curve), 
C2CDist c2c_dist, Curve(*compute_new_curve)(Curves const& curves, Cluster const& cluster));
Points matching_of_vertices(Curve curve_1, Curve curve_2);

bool wedge_parameter_search(Curves const& curves, Clustering& clustering, distance_t(*dist_func)(Curve, Curve), 
C2CDist c2c_dist, distance_t eps, int radius);