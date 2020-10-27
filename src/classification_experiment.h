#ifndef CLASSIFICATION_EXPERIMENT_H
#define CLASSIFICATION_EXPERIMENT_H

#include <functional>
#include <vector>

#include "clustering/center_algs.h"
#include "clustering/init_clustering_algs.h"
#include "Curve.h"
#include "distance_functions.h"

namespace classification {
    using Curves = std::vector<Curve>;

    /**
     * \brief Sample curves with a given period.
     * \param curves The full curve set.
     * \param period The ith curve is selected if i mod period = 0.
     * \return The sampled curves.
     */
    Curves sample(Curves const& curves, unsigned period);

    /**
     * \brief Run the character classification experiment, sampling the curves
     * and using the supplied clustering algorithm.
     *
     * Here we attempt to use clustering for classification: similar shapes
     * should end up in the same cluster and the centre should be a good
     * representative for this class of shapes. Given the handwritten character
     * dataset, we cluster the curves and test whether the clusters correspond
     * to letters. We use 5-fold cross-validation for robustness. We report the
     * results on stdout at the end.
     * \param cluster_alg The basic clustering algorithm.
     * \param center_alg The iterative centre improvement algorithm.
     */
    void classify_characters(
        clustering::ClusterAlg cluster_alg = clustering::ClusterAlg::PAM,
        clustering::CenterAlg center_alg = clustering::CenterAlg::cdba,
        std::function<distance_t(Curve const&, Curve const&)> const& init_dist =
            df::integral_frechet,
        std::function<distance_t(Curve const&, Curve const&)> const& dist =
            df::integral_frechet);

    void identify_pigeons(
        clustering::ClusterAlg cluster_alg = clustering::ClusterAlg::PAM,
        clustering::CenterAlg center_alg = clustering::CenterAlg::cdba,
        std::function<distance_t(Curve const&, Curve const&)> const& init_dist =
            df::integral_frechet,
        std::function<distance_t(Curve const&, Curve const&)> const& dist =
            df::integral_frechet);
}
#endif
