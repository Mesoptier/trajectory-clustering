#ifndef CLUSTERING_ALGS_H
#define CLUSTERING_ALGS_H

#include <functional>
#include <string>
#include <vector>

#include "basic_types.h"
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
        DistanceMatrix<distance_t> const& dist_matrix,
        std::function<distance_t(Curve const&, Curve const&)> const& dist,
        std::function<bool(Curve const&, Curve const&, distance_t)> const& lt_simp,
        bool naive_simplification = false);

    Clustering singleLinkage(Curves const& curves, std::size_t k, std::size_t l,
        SymmetricMatrix const& dist_matrix,
        std::function<bool(Curve const&, Curve const&, distance_t)> const& lt_simp,
        bool naive_simplification);

    Clustering completeLinkage(Curves const& curves, std::size_t k, std::size_t l,
        SymmetricMatrix const& dist_matrix,
        std::function<bool(Curve const&, Curve const&, distance_t)> const& lt_simp,
        bool naive_simplification);

    Clustering gonzalez(Curves const& curves, std::size_t k, std::size_t l,
        CurveSimpMatrix const& dist_matrix,
        std::function<distance_t(Curve const&, Curve const&)> const& dist,
        std::function<bool(Curve const&, Curve const&, distance_t)> const& lt_simp,
        bool naive_simplification);

    Clustering pam_with_simplifications(Curves const& curves,
        std::size_t k, std::size_t l, CurveSimpMatrix const& dist_matrix,
        std::function<bool(Curve const&, Curve const&, distance_t)> const& lt_simp,
        bool naive_simplification);

    /**
     * \brief Update the lists of curve ids in each cluster.
     * \param clustering The clustering with correct centers.
     * \param curves The curves clustered.
     * \param dist The distance function to use.
     */
    void updateClustering(Clustering& clustering, Curves const& curves,
        std::function<distance_t(Curve const&, Curve const&)> const& dist);

    // distance_t kMedianCost(Curves const& curves, Clustering const& clustering,
    //     std::function<distance_t(Curve const&, Curve const&)> const& dist);
}
#endif
