#pragma once

#include <iostream>
#include <map>

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

    // for (int i = 0; i < curves.size(); ++i) {
    //     Curve curve = curves[i];
    //     std::fstream curve_file;
    //     curve_file.open("curves/curve" + std::to_string(i) + ".txt", std::fstream::out | std::fstream::trunc);

    //     for (auto& point: curve.get_points()) {
    //         curve_file << point.x << " " << point.y << " \n";
    //     }
    //     curve_file.close();

    //     script << "\"" << "curves/curve" + std::to_string(i) + ".txt" + "\" with linespoints ls 1 lw 0.5 lt rgb \"" +  + "\" ps 0.01, ";
    // }

    std::vector<std::string> colors = {"red", "yellow", "blue", "yellow", "purple"};

    for (int i = 0; i < clustering.size(); ++i) {
        Cluster cluster = clustering[i];
        std::fstream cluster_center_file;
        cluster_center_file.open("curves/center_curve" + std::to_string(i) + ".txt", std::fstream::out | std::fstream::trunc);
        for (auto& point: cluster.center_curve.get_points()) {
            cluster_center_file << point.x << " " << point.y << " \n";
        }

        for (CurveID curve_id: cluster.curve_ids) {
            Curve curve = curves[curve_id];
            std::fstream curve_file;
            curve_file.open("curves/curve" + std::to_string(curve_id) + ".txt", std::fstream::out | std::fstream::trunc);

            for (auto& point: curve.get_points()) {
                curve_file << point.x << " " << point.y << " \n";
            }
            curve_file.close();

            // script << "\"" << "curves/curve" + std::to_string(curve_id) + ".txt" + "\" with linespoints ls 1 lw 1 lt rgb black ps 0.01, ";
            script << "\"" << "curves/curve" + std::to_string(curve_id) + ".txt" + "\" with linespoints ls 1 lw 1 lt rgb \"" + "black" + "\" ps 0.01, ";
        }

        cluster_center_file.close();
        script << "\"" << "curves/center_curve" + std::to_string(i) + ".txt" + "\" with linespoints ls 1 lw 4 lt rgb \"" + colors[i] +"\", ";

    }

    script.close();
}


void ensemble_experiment() {

    Curves curves = io::read_curves("data/characters/data");
    std::cout << "loaded curves...";

    std::map<char, Curves> curves_by_letter = std::map<char, std::vector<Curve>>();

    for (auto& curve: curves) {
        char letter = curve.name()[21];
        std::cout << letter << "\n";
        auto it = curves_by_letter.find(letter);

        if (it == curves_by_letter.end()) {
            Curves new_vector = {curve};
            curves_by_letter.emplace(letter, new_vector);
        }
        else if (curves_by_letter.at(letter).size() < 8) {
            curves_by_letter.at(letter).push_back(curve);
        }
    }

    std::map<char, std::vector<Curve>>::iterator it;
    std::fstream output_file;
    output_file.open("character_experiment.dat", std::fstream::out | std::fstream::trunc);
    output_file << "character\tfMean\tcdba\tcombination\n";
    for (it = curves_by_letter.begin(); it != curves_by_letter.end(); ++it) {
        std::cout << "comparing on the " << it->first << " curves:\n";
        Clustering fmean_clustering = computeCenterClustering(it->second, 1, 10, ClusterAlg::Gonzalez, CenterAlg::fMean, frechet, average_frechet, "", 1);
        Clustering cdba_clustering = computeCenterClustering(it->second, 1, 10, ClusterAlg::Gonzalez, CenterAlg::cdba, frechet, average_frechet, "", 1);
        Clustering ensemble_clustering = computeCenterClustering(it->second, 1, 10, ClusterAlg::Gonzalez, CenterAlg::ensembleMethod1, frechet, average_frechet, "", 1);

        output_file << fmean_clustering[0].cost << "\t" << cdba_clustering[0].cost << "\t" << ensemble_clustering[0].cost << "\n";

        std::cout << "fmean cost: " << fmean_clustering[0].cost << "\n";
        std::cout << "cdba cost: " << cdba_clustering[0].cost << "\n";
        std::cout << "ensemble cost " << ensemble_clustering[0].cost << "\n";
    }

    output_file.close();
}

void center_update_experiments() {
    Curves curves = read_data();
    // std::cout << "loaded curves...\n";

    // Clustering gonzalez_clustering = computeCenterClustering(curves, 3, 10, ClusterAlg::Gonzalez, CenterAlg::cdba, frechet, integral_frechet, "", 1);
    // std::cout << kMedianCost(curves, gonzalez_clustering, average_frechet);
    // plot(gonzalez_clustering, curves, "plot.txt");

    Clustering gonzalez_clustering = computeCenterClustering(curves, 1, 10, ClusterAlg::Gonzalez, CenterAlg::regression, frechet, average_frechet, "", 1);
    // std::cout << kMedianCost(curves, gonzalez_clustering, average_frechet);
    plot(gonzalez_clustering, curves, "plot.txt");
}

void synthetic_curve_experiment() {
    Curves curves = io::read_curves("synthetic_curves");
    Clustering gonzalez_clustering = computeCenterClustering(curves, 1, 15, ClusterAlg::Gonzalez, CenterAlg::fMean, average_frechet, average_frechet, "", 1);
    plot(gonzalez_clustering, curves, "plot.txt");
}


void preliminary_experiments() {
    using namespace clustering;
    Curves curves = read_data();
    std::cout << curves.size() << "\n";
    // std::vector<> dist_funcs {df::dtw, df::average_frechet, df::frechet};
    std::vector<ClusterAlg> cluster_algs = {ClusterAlg::CompleteLinkage, ClusterAlg::Gonzalez, ClusterAlg::PAM};

    Clustering gonzalez_f_mean = computeCenterClustering(curves, 1, 10, ClusterAlg::Gonzalez, CenterAlg::fMean, average_frechet, average_frechet, "", 1);
    std::cout << "computed average frechet clustering\n";
    // Clustering gonzalez_f_center = computeCenterClustering(curves, 4, 10, ClusterAlg::Gonzalez, CenterAlg::fCenter, frechet, "", 1);
    // std::cout << "computed gonzalez frechet clustering\n";
    // Clustering gonzalez_dtw = computeCenterClustering(curves, 4, 10, ClusterAlg::Gonzalez, CenterAlg::dtwMean, dtw, "", 1);
    // std::cout << "computed gonzalez dtw clustering\n";

    distance_t av_frechet_median_cost = 0;
    distance_t av_frechet_center_cost = 0;
    distance_t av_frechet_dtw_median_cost = 0;
    for (auto& cluster: gonzalez_f_mean) {
        av_frechet_median_cost += calcC2CDist(curves, cluster.center_curve, cluster.curve_ids, C2CDist::Median, average_frechet);
        av_frechet_center_cost += calcC2CDist(curves, cluster.center_curve, cluster.curve_ids, C2CDist::Max, frechet);
        av_frechet_dtw_median_cost += calcC2CDist(curves, cluster.center_curve, cluster.curve_ids, C2CDist::Median, dtw);
    }

    // distance_t frechet_median_cost = 0;
    // distance_t frechet_center_cost = 0;
    // distance_t frechet_dtw_median_cost = 0;
    // for (auto& cluster: gonzalez_f_center) {
    //     frechet_median_cost += calcC2CDist(curves, cluster.center_curve, cluster.curve_ids, C2CDist::Median, average_frechet);
    //     frechet_center_cost += calcC2CDist(curves, cluster.center_curve, cluster.curve_ids, C2CDist::Max, frechet);
    //     frechet_dtw_median_cost += calcC2CDist(curves, cluster.center_curve, cluster.curve_ids, C2CDist::Median, dtw);
    // }

    // distance_t dtw_av_frechet_median_cost = 0;
    // distance_t dtw_frechet_center_cost = 0;
    // distance_t dtw_median_cost = 0;
    // for (auto& cluster: gonzalez_dtw) {
    //     dtw_av_frechet_median_cost += calcC2CDist(curves, cluster.center_curve, cluster.curve_ids, C2CDist::Median, average_frechet);
    //     dtw_frechet_center_cost += calcC2CDist(curves, cluster.center_curve, cluster.curve_ids, C2CDist::Max, frechet);
    //     dtw_median_cost += calcC2CDist(curves, cluster.center_curve, cluster.curve_ids, C2CDist::Median, dtw);
    // }

    std::cout << "gonzalez + average frechet mean results:\n";
    std::cout << "av_frechet_median_cost:\t" << av_frechet_median_cost << "\n";
    std::cout << "av_frechet_center_cost:\t" << av_frechet_center_cost << "\n";
    std::cout << "av_frechet_dtw_median_cost:\t" << av_frechet_dtw_median_cost << "\n\n";

    // std::cout << "gonzalez + frechet centering results:\n";
    // std::cout << "frechet_median_cost:\t" << frechet_median_cost << "\n";
    // std::cout << "frechet_center_cost:\t" << frechet_center_cost << "\n";
    // std::cout << "frechet_dtw_median_cost:\t" << frechet_dtw_median_cost << "\n\n";

    // std::cout << "gonzalez + dtw results:\n";
    // std::cout << "dtw_av_frechet_median_cost:\t" << dtw_av_frechet_median_cost << "\n";
    // std::cout << "dtw_frechet_center_cost:\t" << dtw_frechet_center_cost << "\n";
    // std::cout << "dtw_median_cost:\t" << dtw_median_cost << "\n\n";


    // Clustering pam = computeCenterClustering(curves, 4, 10, ClusterAlg::Pam, CenterAlg::fMean, average_frechet, "pigeon_matrix.txt", 1);
    // Clustering gonzalez_frechet = computeCenterClustering(curves, 4, 10, ClusterAlg::Gonzalez, CenterAlg::fCenter, frechet, "", 1);
    // Clustering gonzalez_mean = computeCenterClustering(curves, 4, 10, ClusterAlg::Gonzalez, CenterAlg::fMean, average_frechet, "", 1);
    // Clustering gonzalez_average_center = computeCenterClustering(curves, 4, 10, ClusterAlg::Gonzalez, CenterAlg::fCenter, average_frechet, "", 1);
    // Clustering pam_dtw = computeCenterClustering(curves, 4, 10, ClusterAlg::Pam, CenterAlg::dtwMean, dtw, "pigeon_matrix_dtw.txt", 1);
    // Clustering gonzalez_mean_dtw = computeCenterClustering(curves, 4, 10, ClusterAlg::Gonzalez, CenterAlg::dtwMean, dtw, "", 1);

    // std::cout << "pam + frechet mean: " << kMedianCost(curves, pam, average_frechet) << "\n";
    // std::cout << "gonzalez frechet: " << kMedianCost(curves, gonzalez_frechet, average_frechet) << "\n";
    // std::cout << "gonzalez mean: " << kMedianCost(curves, gonzalez_mean, average_frechet) << "\n";
    // std::cout << "gonzalez average center: " << kMedianCost(curves, gonzalez_average_center, average_frechet) << "\n";
    // std::cout << "pam + dtw mean centering: " << kMedianCost(curves, pam_dtw, average_frechet) << "\n";
    // std::cout << "gonzalez mean dtw: " << kMedianCost(curves, gonzalez_mean_dtw, average_frechet) << "\n"; 
}

void wedge_method_experiment() {
    Curves curves = io::read_curves("data/characters/data");
    std::cout << "loaded curves...";

    std::map<char, Curves> curves_by_letter = std::map<char, std::vector<Curve>>();

    for (auto& curve: curves) {
        char letter = curve.name()[21];
        std::cout << letter << "\n";
        auto it = curves_by_letter.find(letter);

        if (it == curves_by_letter.end()) {
            Curves new_vector = {curve};
            curves_by_letter.emplace(letter, new_vector);
        }
        else if (curves_by_letter.at(letter).size() < 11) {
            curves_by_letter.at(letter).push_back(curve);
        }
    }

    std::map<char, std::vector<Curve>>::iterator it;

    std::fstream output_file;
    output_file.open("wedge_experiment.dat", std::fstream::out | std::fstream::trunc);

    output_file << "character\tcdba_cost\twedge_cost\n";

    for (it = curves_by_letter.begin(); it != curves_by_letter.end(); ++it) {
        std::cout << "comparing the " << it->first << " curves...\n";
        Clustering cdba_clustering = computeCenterClustering(it->second, 1, 12, ClusterAlg::Gonzalez, CenterAlg::cdba, frechet, average_frechet, "", 1);
        Clustering wedge_clustering = computeCenterClustering(it->second, 1, 12, ClusterAlg::Gonzalez, CenterAlg::wedge, frechet, average_frechet, "", 1);
        output_file << it->first << "\t" << std::to_string(cdba_clustering[0].cost) << "\t" << std::to_string(wedge_clustering[0].cost) << "\n";
    }

    output_file.close();
}
}
