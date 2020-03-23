#pragma once

#include "Curve.h"
#include "distance_functions.h"
#include "clustering/clustering_algs.h"
#include "clustering/center_algs.h"
#include "clustering/center_clustering_algs.h"
#include "io.h"

Curves read_data() {
    Curves curves = io::read_pigeon_curves("data/Data_for_Mann_et_al_RSBL 2/Bladon & Church route recapping/bladon heath");
    Curves simplified_curves = Curves();
    for (auto curve: curves) {
        simplified_curves.push_back(curve.naive_l_simplification(50));
    }
    return simplified_curves;
}

void preliminary_experiments() {
    Curves curves = read_data();
    std::cout << curves.size() << "\n";
    std::vector<distance_t (*)(Curve, Curve)> dist_funcs = {dtw, average_frechet, frechet};
    std::vector<ClusterAlg> cluster_algs = {
		ClusterAlg::CompleteLinkage, ClusterAlg::Gonzalez, ClusterAlg::Pam};

    Clustering gonzalez_frechet = computeCenterClustering(curves, 5, 10, ClusterAlg::Gonzalez, CenterAlg::fCenter, frechet, 1);
    std::cout << gonzalez_frechet.size() << "\n";
}