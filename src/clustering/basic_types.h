#ifndef BASIC_TYPES_H
#define BASIC_TYPES_H

#include <limits>
#include <vector>

#include "Curve.h"
#include "geom.h"
#include "utils/id.h"

struct Cluster {
    CurveIDs curve_ids;
    Curve center_curve;
    CurveID center_id;
    distance_t cost = std::numeric_limits<distance_t>::max();
};

using Clustering = std::vector<Cluster>;
using ClusterID = ID<Cluster>;
using ClusterIDs = std::vector<ClusterID>;
#endif
