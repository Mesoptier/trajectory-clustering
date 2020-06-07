#pragma once

#include <iostream>

#include "Curve.h"
#include "distance_functions.h"
#include "clustering/clustering_algs.h"
#include "clustering/center_algs.h"
#include "clustering/center_clustering_algs.h"
#include "utils/io.h"

namespace experiments {
Curves read_data() {
    // Curves curves = io::read_pigeon_curves("data/Data_for_Mann_et_al_RSBL 2/Bladon & Church route recapping/bladon heath");
    // Curves simplified_curves = Curves();
    // for (auto curve: curves) {
    //     simplified_curves.push_back(curve.naive_l_simplification(50));
    // }
    // return simplified_curves;

    Curves curves = io::read_curves("data/characters/data");
    Curves ret(curves.begin(), curves.begin() + 70);
    // Curves bs = Curves(curves.begin() + 200, curves.begin() + 210);

    // Curves final_curves = Curves();
    // for (auto& curve: as) {
    //     final_curves.push_back(curve);
    // }
    // for (auto& curve: bs) {
        // final_curves.push_back(curve);
    // }
    return ret;
}

void plot(Clustering& clustering, Curves& curves, std::string script_name) {
    std::fstream script;
    script.open(script_name, std::fstream::out | std::fstream::trunc);
    script << "plot ";
    for (std::size_t i = 0; i < curves.size(); ++i) {
        Curve curve = curves[i];
        std::fstream curve_file;
        curve_file.open("curves/curve_" + std::to_string(i) + ".txt", std::fstream::out | std::fstream::trunc);

        for (auto& point: curve.get_points()) {
            curve_file << point.x << " " << point.y << " \n";
        }
        curve_file.close();

        script << "\"" << "curves/curve_" + std::to_string(i) + ".txt" + "\" with linespoints ls 1 lw 0.5 lt rgb \"black\" ps 0.01, ";
    }

    std::vector<std::string> colors = {"red", "green", "blue", "yellow", "purple"};

    for (std::size_t i = 0; i < clustering.size(); ++i) {
        Cluster cluster = clustering[i];
        std::fstream cluster_center_file;
        cluster_center_file.open("curves/center_" + std::to_string(i) + ".txt", std::fstream::out | std::fstream::trunc);
        for (auto& point: cluster.center_curve.get_points()) {
            cluster_center_file << point.x << " " << point.y << " \n";
        }

        cluster_center_file.close();
        script << "\"" << "curves/center_" + std::to_string(i) + ".txt" + "\" with linespoints ls 1 lw 2 lt rgb \"" + colors[i] +"\", ";
    }

    script.close();
    

}

void center_update_experiments() {
    Curves curves = read_data();
    // std::cout << "loaded curves...\n";

    Clustering gonzalez_clustering = clustering::computeCenterClustering(curves, 1, 15, clustering::ClusterAlg::Gonzalez, clustering::CenterAlg::ensembleMethod1, df::average_frechet, "", 1);
    // std::cout << kMedianCost(curves, gonzalez_clustering, average_frechet);
    plot(gonzalez_clustering, curves, "test_script.txt");
}

void synthetic_curve_experiment() {
    Curves curves = io::read_curves("synthetic_curves");
    Clustering gonzalez_clustering = clustering::computeCenterClustering(curves, 1, 15, clustering::ClusterAlg::Gonzalez, clustering::CenterAlg::fMean, df::average_frechet, "", 1);
    plot(gonzalez_clustering, curves, "test_script.txt");
}


void preliminary_experiments() {
    using namespace clustering;
    Curves curves = read_data();
    std::cout << curves.size() << "\n";
    // std::vector<> dist_funcs {df::dtw, df::average_frechet, df::frechet};
    std::vector<ClusterAlg> cluster_algs = {ClusterAlg::CompleteLinkage, ClusterAlg::Gonzalez, ClusterAlg::PAM};

    Clustering pam = computeCenterClustering(curves, 4, 10, ClusterAlg::PAM, CenterAlg::fMean, df::average_frechet, "pigeon_matrix.txt", 1);
    Clustering gonzalez_frechet = computeCenterClustering(curves, 4, 10, ClusterAlg::Gonzalez, CenterAlg::fCenter, df::frechet, "", 1);
    Clustering gonzalez_mean = computeCenterClustering(curves, 4, 10, ClusterAlg::Gonzalez, CenterAlg::fMean, df::average_frechet, "", 1);
    Clustering gonzalez_average_center = computeCenterClustering(curves, 4, 10, ClusterAlg::Gonzalez, CenterAlg::fCenter, df::average_frechet, "", 1);
    Clustering pam_dtw = computeCenterClustering(curves, 4, 10, ClusterAlg::PAM, CenterAlg::dtwMean, df::dtw, "pigeon_matrix_dtw.txt", 1);
    Clustering gonzalez_mean_dtw = computeCenterClustering(curves, 4, 10, ClusterAlg::Gonzalez, CenterAlg::dtwMean, df::dtw, "", 1);

    std::cout << "pam + frechet mean: " << kMedianCost(curves, pam, df::average_frechet) << "\n";
    std::cout << "gonzalez frechet: " << kMedianCost(curves, gonzalez_frechet, df::average_frechet) << "\n";
    std::cout << "gonzalez mean: " << kMedianCost(curves, gonzalez_mean, df::average_frechet) << "\n";
    std::cout << "gonzalez average center: " << kMedianCost(curves, gonzalez_average_center, df::average_frechet) << "\n";
    std::cout << "pam + dtw mean centering: " << kMedianCost(curves, pam_dtw, df::average_frechet) << "\n";
    std::cout << "gonzalez mean dtw: " << kMedianCost(curves, gonzalez_mean_dtw, df::average_frechet) << "\n"; 
}
}
