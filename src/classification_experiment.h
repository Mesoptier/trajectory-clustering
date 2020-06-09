#ifndef CLASSIFICATION_EXPERIMENT_H
#define CLASSIFICATION_EXPERIMENT_H

#include <vector>

#include "clustering/center_algs.h"
#include "clustering/clustering_algs.h"
#include "Curve.h"

namespace classification {
    using Curves = std::vector<Curve>;

    /**
     * \brief Sample curves with a given period.
     * \param curves The full curve set.
     * \param period The ith curve is selected if i mod period = 0.
     * \return The sampled curves.
     */
    Curves sample(const Curves& curves, unsigned period);

    /**
     * \brief Run the character classification experiment, sampling the curves
     * and using the supplied clustering algorithm.
     *
     * Here we attempt to use clustering for classification: similar shapes
     * should end up in the same cluster and the centre should be a good
     * representative for this class of shapes. Given the handwritten character
     * dataset, we cluster the curves and test whether the clusters correspond
     * to letters. We use 5-fold cross-validation for robustness. We report the
     * confusion matrix on stdout at the end.
     * \param cluster_alg The basic clustering algorithm.
     * \param center_alg The iterative centre improvement algorithm.
     */
    void characterClassification(
        clustering::ClusterAlg cluster_alg = clustering::ClusterAlg::Gonzalez,
        clustering::CenterAlg center_alg = clustering::CenterAlg::fCenter);
}
#endif
