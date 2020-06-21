#ifndef BASIC_TYPES_H
#define BASIC_TYPES_H

#include <vector>

#include "Curve.h"
#include "geom.h"
#include "id.h"

struct Cluster
{
    CurveIDs curve_ids;
    Curve center_curve;
    distance_t cost = std::numeric_limits<distance_t>::max();
    CurveID center_id;
};

using Clustering = std::vector<Cluster>;
using ClusterID = ID<Cluster>;
using ClusterIDs = std::vector<ClusterID>;
#endif
