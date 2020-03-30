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

    // Curves curves = io::read_curves("data/characters/data");
    // return Curves(curves.begin(), curves.begin() + 30);
}


void preliminary_experiments() {
    Curves curves = read_data();
    std::cout << curves.size() << "\n";
    std::vector<distance_t (*)(Curve, Curve)> dist_funcs = {dtw, average_frechet, frechet};
    std::vector<ClusterAlg> cluster_algs = {ClusterAlg::CompleteLinkage, ClusterAlg::Gonzalez, ClusterAlg::Pam};

    Clustering pam = computeCenterClustering(curves, 4, 10, ClusterAlg::Pam, CenterAlg::fMean, average_frechet, "pigeon_matrix.txt", 1);
    Clustering gonzalez_frechet = computeCenterClustering(curves, 4, 10, ClusterAlg::Gonzalez, CenterAlg::fCenter, frechet, "", 1);
    Clustering gonzalez_mean = computeCenterClustering(curves, 4, 10, ClusterAlg::Gonzalez, CenterAlg::fMean, average_frechet, "", 1);
    Clustering gonzalez_average_center = computeCenterClustering(curves, 4, 10, ClusterAlg::Gonzalez, CenterAlg::fCenter, average_frechet, "", 1);
    Clustering pam_dtw = computeCenterClustering(curves, 4, 10, ClusterAlg::Pam, CenterAlg::dtwMean, dtw, "pigeon_matrix_dtw.txt", 1);
    Clustering gonzalez_mean_dtw = computeCenterClustering(curves, 4, 10, ClusterAlg::Gonzalez, CenterAlg::dtwMean, dtw, "", 1);

    std::cout << "pam + frechet mean: " << kMedianCost(curves, pam, average_frechet) << "\n";
    std::cout << "gonzalez frechet: " << kMedianCost(curves, gonzalez_frechet, average_frechet) << "\n";
    std::cout << "gonzalez mean: " << kMedianCost(curves, gonzalez_mean, average_frechet) << "\n";
    std::cout << "gonzalez average center: " << kMedianCost(curves, gonzalez_average_center, average_frechet) << "\n";
    std::cout << "pam + dtw mean centering: " << kMedianCost(curves, pam_dtw, average_frechet) << "\n";
    std::cout << "gonzalez mean dtw: " << kMedianCost(curves, gonzalez_mean_dtw, average_frechet) << "\n"; 
}