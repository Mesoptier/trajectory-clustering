#ifndef CLASSIFICATION_EXPERIMENT_H
#define CLASSIFICATION_EXPERIMENT_H

#include <vector>

#include "clustering/center_algs.h"
#include "clustering/clustering_algs.h"
#include "Curve.h"

namespace classification {
    using Curves = std::vector<Curve>;

    Curves sample(const Curves& curves, unsigned period);

    void characterClassification(
        clustering::ClusterAlg cluster_alg = clustering::ClusterAlg::Gonzalez,
        clustering::CenterAlg center_alg = clustering::CenterAlg::fCenter);
}
#endif
