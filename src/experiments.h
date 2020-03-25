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

void plot(Clustering& clustering, Curves& curves, std::string script_name) {
    std::fstream script;
    script.open(script_name, std::fstream::out | std::fstream::trunc);
    script << "plot ";
    for (int i = 0; i < curves.size(); ++i) {
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

    for (int i = 0; i < clustering.size(); ++i) {
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

void preliminary_experiments() {
    Curves curves = read_data();
    std::cout << curves.size() << "\n";
    std::vector<distance_t (*)(Curve, Curve)> dist_funcs = {dtw, average_frechet, frechet};
    std::vector<ClusterAlg> cluster_algs = {
		ClusterAlg::CompleteLinkage, ClusterAlg::Gonzalez, ClusterAlg::Pam};

    Clustering gonzalez_frechet = computeCenterClustering(curves, 5, 10, ClusterAlg::Gonzalez, CenterAlg::fCenter, frechet, 1);
}