#ifndef CLUSTERING_ALGS_H
#define CLUSTERING_ALGS_H

#include <functional>
#include <string>
#include <vector>

#include "basic_types.h"
#include "clustering/center_algs.h"
#include "clustering/pam.h"
#include "utils/CurveSimpMatrix.h"

namespace clustering {
    using Curves = std::vector<Curve>;

    enum class ClusterAlg {
        SingleLinkage,
        CompleteLinkage,
        Gonzalez,
        PAM
    };
    std::string toString(ClusterAlg cluster_alg);

    Clustering computeClustering(Curves const& curves,
        std::size_t k, std::size_t l, ClusterAlg cluster_alg,
        std::function<distance_t(Curve const&, Curve const&)> const& dist,
        std::string const& dist_matrix, bool naive_simplification = false);

    Clustering singleLinkage(Curves const& curves,
        std::size_t k, std::size_t l,
        std::function<distance_t(Curve const&, Curve const&)> const& dist,
        bool naive_simplification);

    Clustering completeLinkage(Curves const& curves,
        std::size_t k, std::size_t l,
        std::function<distance_t(Curve const&, Curve const&)> const& dist,
        bool naive_simplification);

    Clustering gonzalez(Curves const& curves,
        std::size_t k, std::size_t l,
        std::function<distance_t(Curve const&, Curve const&)> const& dist,
        bool naive_simplification);

    // assign curves to closest clusters
    void updateClustering(Curves const& curves, Clustering& clustering,
        std::function<distance_t(Curve const&, Curve const&)> const& dist);

    distance_t calcDiameter(Curves const& curves,
        std::function<distance_t(Curve const&, Curve const&)> const& dist);

    distance_t kMedianCost(Curves const& curves, Clustering const& clustering,
        std::function<distance_t(Curve const&, Curve const&)> const& dist);

    Clustering pam_with_simplifications(Curves const& curves,
        std::size_t k, std::size_t l,
        std::function<distance_t(Curve const&, Curve const&)> const& dist,
        std::string const& matrix_file_name = "",
        bool naive_simplification = false);

    Clustering pam_with_centering(Curves const& curves,
        std::size_t k, std::size_t l,
        std::function<distance_t(Curve const&, Curve const&)> const& dist,
        std::string const& matrix_file_name = "",
        bool naive_simplification = false);
}
#endif
