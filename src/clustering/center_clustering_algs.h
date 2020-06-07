#ifndef CENTER_CLUSTERING_ALGS_H
#define CENTER_CLUSTERING_ALGS_H

#include <cstddef>
#include <functional>
#include <string>
#include "clustering/center_algs.h"
#include "clustering/clustering_algs.h"
#include "Curve.h"

namespace clustering {
    Clustering computeCenterClustering(Curves const& curves,
        std::size_t k, std::size_t l,
        ClusterAlg cluster_alg, CenterAlg center_alg,
        std::function<distance_t(Curve const&, Curve const&)> const& dist_func,
        std::string const& dist_matrix,
        unsigned max_rounds = 10);
}
#endif
