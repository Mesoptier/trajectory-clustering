#pragma once

#include "Curve.h"
#include "distance_functions.h"
#include "clustering/clustering_algs.h"
#include "clustering/center_algs.h"
#include "clustering/center_clustering_algs.h"
#include "io.h"
#include "synthetic_curves.h"
#include "simplification/imaiiri.h"
#include <map>
#include "clustering/plot_clustering.h"
// #include "ortools/linear_solver/linear_solver.h"


Curves read_data() {
    // Curves curves = io::read_pigeon_curves("data/Data_for_Mann_et_al_RSBL 2/Bladon & Church route recapping/bladon heath/c35");
    // Curves simplified_curves = Curves();
    // // for (auto curve: curves) {
    // //     // std::cout << curve.size() << "\n";
    // //     simplified_curves.push_back(curve.naive_l_simplification(50));
    // // }
    // // return simplified_curves;
    // return curves;

    Curves curves = io::read_curves("data/characters/data");
    Curves as = Curves(curves.begin(), curves.begin() + 5);
    return as;
    // Curves bs = Curves(curves.begin() + 200, curves.begin() + 210);

    // Curves final_curves = Curves();
    // for (auto& curve: as) {
    //     final_curves.push_back(curve);
    // }
    // // for (auto& curve: bs) {
    //     // final_curves.push_back(curve);
    // // }
    // return final_curves;
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

void synthetic_curve_experiment() {
    Curves curves = io::read_curves("synthetic_curves");
    Clustering gonzalez_clustering = computeCenterClustering(curves, 1, 15, ClusterAlg::Gonzalez, CenterAlg::fMean, average_frechet, average_frechet, "", 1);
    plot_clustering(gonzalez_clustering, curves, "plot.txt");
}


void preliminary_experiments() {

    Curves curves = read_data();
    std::cout << curves.size() << "\n";
    std::vector<distance_t (*)(Curve, Curve)> dist_funcs = {dtw, average_frechet, frechet};
    std::vector<ClusterAlg> cluster_algs = {ClusterAlg::CompleteLinkage, ClusterAlg::Gonzalez, ClusterAlg::Pam};


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
        else if (curves_by_letter.at(letter).size() < 21) {
            curves_by_letter.at(letter).push_back(curve);
        }
    }

    std::map<char, std::vector<Curve>>::iterator it;

    std::fstream output_file;
    output_file.open("results/cdba_wedge_character_experiment/wedge_experiment.dat", std::fstream::out | std::fstream::trunc);

    std::fstream plot_file;
    plot_file.open("results/cdba_wedge_character_experiment/plots.txt");
    // output_file << "character\tcdba_cost\twedge_cost\n";

    for (it = curves_by_letter.begin(); it != curves_by_letter.end(); ++it) {
        std::cout << "comparing the " << it->first << " curves...\n";
        Clustering cdba_clustering = computeCenterClustering(it->second, 1, 10, ClusterAlg::Gonzalez, CenterAlg::cdba, frechet, integral_frechet<5>, "", 3);
        Clustering wedge_clustering = computeCenterClustering(it->second, 1, 10, ClusterAlg::Gonzalez, CenterAlg::wedge, frechet, integral_frechet<5>, "", 3);
        std::cout << "cdba_cost" << cdba_clustering[0].cost << "\n";
        std::cout << "wedge clustering" << wedge_clustering[0].cost << "\n";
        plot_clustering(cdba_clustering, it->second, "results/cdba_wedge_character_experiment/cdba_" + std::string(1, it->first) + ".txt");
        plot_clustering(wedge_clustering, it->second, "results/cdba_wedge_character_experiment/wedge_" + std::string(1, it->first) + ".txt");
        output_file << it->first << "\t" << std::to_string(cdba_clustering[0].cost) << "\t" << std::to_string(wedge_clustering[0].cost) << "\n";
        plot_file << "cdba_" + std::string(1, it->first) << std::endl;
        plot_file << "wedge_" + std::string(1, it->first) << std::endl;
    }
    output_file.close();
    plot_file.close();
}

void pigeon_experiment() {
    std::vector<std::string> pigeons = {
        "a55", "brc", "c17", "c35",
        "p29", "p39", "p94"
    };

    std::fstream clustering_names;
    clustering_names.open("results/pigeon_visual_experiment/plots.txt", std::fstream::out | std::fstream::trunc);

    for (std::string pigeon: pigeons) {
        Curves raw_curves = io::read_pigeon_curves_utm("data/Data_for_Mann_et_al_RSBL 2/Bladon & Church route recapping/bladon heath/" + pigeon);
        Curves curves = Curves();
        for (auto curve: raw_curves) {
            curves.push_back(curve.naive_l_simplification(100));
        }


        std::cout << "dba...\n";
        Clustering dba_clustering = computeCenterClustering(curves, 3, 10, ClusterAlg::Gonzalez, CenterAlg::dba, frechet, dtw, "", 1);
        std::cout << "fsa...\n";
        Clustering fsa_clustering = computeCenterClustering(curves, 3, 10, ClusterAlg::Gonzalez, CenterAlg::fCenter, frechet, frechet, "", 1);
        // std::cout << "cdba av...\n";
        // Clustering cdba_average_clustering = computeCenterClustering(curves, 1, 10, ClusterAlg::Gonzalez, CenterAlg::cdba, frechet, average_frechet, "", 5);
        std::cout << "cdba int...\n";
        Clustering cdba_integral_clustering = computeCenterClustering(curves, 3, 10, ClusterAlg::Gonzalez, CenterAlg::cdba, frechet, integral_frechet<250>, "", 1);
        // std::cout << "wedge av...\n";
        // Clustering wedge_average_clustering = computeCenterClustering(curves, 1, 10, ClusterAlg::Gonzalez, CenterAlg::wedge, frechet, average_frechet, "", 5);
        std::cout << "wedge int...\n";
        Clustering wedge_integral_clustering = computeCenterClustering(curves, 3, 10, ClusterAlg::Gonzalez, CenterAlg::wedge, frechet, integral_frechet<250>, "", 1);

        plot_clustering(dba_clustering, curves, "results/pigeon_visual_experiment/" + pigeon + "_dba.txt");
        plot_clustering(fsa_clustering, curves, "results/pigeon_visual_experiment/" + pigeon + "_fsa.txt");
        // plot_clustering(cdba_average_clustering, curves, "results/pigeon_visual_experiment/" + pigeon + "_cdba_average_1.txt");
        plot_clustering(cdba_integral_clustering, curves, "results/pigeon_visual_experiment/" + pigeon + "_cdba_integral.txt");
        // plot_clustering(wedge_average_clustering, curves, "results/pigeon_visual_experiment/" + pigeon + "_wedge_average.txt");
        plot_clustering(wedge_integral_clustering, curves, "results/pigeon_visual_experiment/" + pigeon + "_wedge_integral.txt");

        clustering_names << pigeon + "_dba\n";
        clustering_names << pigeon + "_fsa\n";
        // clustering_names << pigeon + "_cdba_average_1\n";
        clustering_names << pigeon + "_cdba_integral\n";
        // clustering_names << pigeon + "_wedge_average\n";
        clustering_names << pigeon + "_wedge_integral\n";
    }

    clustering_names.close();
}



void find_params() {
    Curves raw_curves = io::read_pigeon_curves_utm("data/Data_for_Mann_et_al_RSBL 2/Bladon & Church route recapping/bladon heath/a55");
    std::cout << raw_curves[0][0];
    
    Curves curves = Curves();
    for (auto curve: raw_curves) {
        curves.push_back(curve.naive_l_simplification(100));
    }

    Clustering initial_clustering = computeClustering(curves, 3, 10, ClusterAlg::Gonzalez, frechet, "", false);

    std::vector<distance_t> epsilons = {10, 30, 50, 70};
    std::vector<int> radii = {5, 10, 20, 30, 40};

    for (distance_t eps: epsilons) {
        for (int radius: radii) {
            std::cout << "attempting with eps=" + std::to_string(eps) + " and radius=" + std::to_string(radius) << "\n";
            wedge_parameter_search(curves, initial_clustering, integral_frechet<5>, C2CDist::Median, eps, radius);
        }
    }
}

void full_clustering_experiment() {
    std::vector<std::string> pigeon_sites = {
        "Bladon & Church route recapping/bladon heath",
        "Bladon & Church route recapping/church hanborough",
        "Horspath",
        "Weston"
    };

    std::vector<std::string> dist_matrix_paths = {
        "bladon_heath_10.txt",
        "church_10.txt",
        "horspath_10.txt",
        "weston_10.txt"
    };

    std::fstream results;
    results.open("results/initial_clustering_experiment/results.txt", std::fstream::out | std::fstream::trunc);

    for (int i = 0 ; i < pigeon_sites.size(); ++i) {
        std::string site = pigeon_sites[i];
        std::string dist_matrix_path = dist_matrix_paths[i];
        CurveSimpMatrix dist_matrix = CurveSimpMatrix(dist_matrix_path);

        Curves raw_pigeon_curves = io::read_pigeon_curves_utm("data/Data_for_Mann_et_al_RSBL 2/" + site);
        Curves curves = Curves();
        Curves simplifications = Curves();

        int count = 0;
        for (auto curve: raw_pigeon_curves) {
            Curve naive_simp = curve.naive_l_simplification(100);
            curves.push_back(naive_simp);
            // simplifications.push_back(simplify(naive_simp, 10, frechet));
            // ++count;
            // std::cout << count << std::endl;
        }

        for (int k = 1; k <= 10; ++k) {
            std::cout << site << "\n";
            Clustering gonzalez_pam_clustering = computeClustering(
                curves, k, 10, ClusterAlg::GonzalezPam, integral_frechet<250>,
                dist_matrix_path, false
            );
            for (int i = 0; i < 5; ++i) {
                Clustering alt = computeClustering(
                    curves, k, 10, ClusterAlg::GonzalezPam, integral_frechet<250>,
                    dist_matrix_path, false
                );

                if (kMedianCostMat(curves, alt, dist_matrix) < kMedianCostMat(curves, gonzalez_pam_clustering, dist_matrix)) {
                    gonzalez_pam_clustering = std::move(alt);
                }
            }

            Clustering pam_clustering = computeClustering(
                curves, k, 10, ClusterAlg::Pam, integral_frechet<250>,
                dist_matrix_path, false
            );

            Clustering gonzalez_clustering = computeClustering(
                curves, k, 10, ClusterAlg::Gonzalez, integral_frechet<250>,
                dist_matrix_path, false
            );
            for (int i = 0; i < 5; ++i) {
                Clustering alt = computeClustering(
                    curves, k, 10, ClusterAlg::Gonzalez, integral_frechet<250>,
                    dist_matrix_path, false
                );

                if (kMedianCostMat(curves, alt, dist_matrix) < kMedianCostMat(curves, gonzalez_clustering, dist_matrix)) {
                    gonzalez_clustering = std::move(alt);
                }
            }


            results << site << "\t" << k << "\t" <<
            kMedianCostMat(curves, gonzalez_pam_clustering, dist_matrix) <<
            "\t" << kMedianCostMat(curves, pam_clustering, dist_matrix) <<
            "\t" << kMedianCostMat(curves, gonzalez_clustering, dist_matrix) << std::endl;
        }
    }

    results.close();
}

void compute_matrices() {

    Curves raw_pigeon_curves = io::read_pigeon_curves_utm("data/Data_for_Mann_et_al_RSBL 2/Weston");
    // Curves raw_characters = io::read_curves("data/characters/data");

    Curves curves = Curves();
    Curves simplifications = Curves();

    int count = 0;
    for (auto curve: raw_pigeon_curves) {
        std::cout << curve.name() << std::endl;
        Curve naive_simp = curve.naive_l_simplification(100);
        curves.push_back(naive_simp);
        simplifications.push_back(simplify(naive_simp, 10, frechet));
        ++count;
        std::cout << count << std::endl;
    }

    CurveSimpMatrix pigeon_matrix = CurveSimpMatrix(curves, simplifications, integral_frechet<250>);
    pigeon_matrix.write("weston_10.txt");
}

distance_t kMedianCost(Curves const& curves, Clustering const& clustering, distance_t(*dist_func)(Curve, Curve)) {
	distance_t cost = 0;

	for (auto& cluster: clustering) {
		cost += calcC2CDist(curves, cluster.center_curve, cluster.curve_ids, C2CDist::Median, dist_func);
	}

	return cost;
}

void center_update_experiment_characters(std::string directory, int n, int k, int l) {
    Curves curves = io::read_curves("data/characters/data");
    std::cout << "loaded curves...";

    std::map<char, Curves> curves_by_letter = std::map<char, std::vector<Curve>>();

    for (auto& curve: curves) {
        char letter = curve.name()[21];
        auto it = curves_by_letter.find(letter);

        if (it == curves_by_letter.end()) {
            Curves new_vector = {curve};
            curves_by_letter.emplace(letter, new_vector);
        }
        else if (curves_by_letter.at(letter).size() < n) {
            curves_by_letter.at(letter).push_back(curve);
        }
    }

    std::map<char, std::vector<Curve>>::iterator it;

    system(("mkdir results/" + directory).c_str());

    std::fstream plots;
    plots.open("results/" + directory + "/plots.txt", std::fstream::out | std::fstream::trunc);

    std::fstream k_med_scores;
    k_med_scores.open("results/" + directory + "/k_med_scores.txt", std::fstream::out | std::fstream::trunc);

    k_med_scores << "characters,DBA,FSA,CDBA,WEDGE\n";

    for (it = curves_by_letter.begin(); it != curves_by_letter.end(); ++it) {
        std::cout << it->first << "...\n";

        Curves curves = it->second;

        Clustering initial_clustering = computeClustering(curves, k, l, ClusterAlg::Pam, integral_frechet<5>, "", false);

        Clustering dba_clustering = initial_clustering;
        int count = 1;
        int max_count = 20;
        while (count <= max_count && computerCenters(curves, dba_clustering, l, CenterAlg::dbaChar, dtw)) {
            updateClustering(curves, dba_clustering, dtw, nullptr);
            ++count;
        }

        Clustering fsa_clustering = initial_clustering;
        count = 1;
        while (count <= max_count && computerCenters(curves, fsa_clustering, l, CenterAlg::fCenter, frechet)) {
            updateClustering(curves, fsa_clustering, frechet, nullptr);
            ++count;
        }

        Clustering cdba_clustering = initial_clustering;
        count = 1;
        while (count <= max_count && computerCenters(curves, cdba_clustering, l, CenterAlg::cdbaChar, integral_frechet<5>)) {
            updateClustering(curves, cdba_clustering, integral_frechet<5>, nullptr);
            ++count;
        }

        Clustering wedge_clustering = initial_clustering;
        count = 1;
        while (count <= max_count && computerCenters(curves, wedge_clustering, l, CenterAlg::wedgeChar, integral_frechet<5>)) {
            updateClustering(curves, wedge_clustering, integral_frechet<5>, nullptr);
            ++count;
        }

        plot_clustering(dba_clustering, it->second, "results/" + directory + "/dba_" + std::string(1, it->first) + ".txt");
        plot_clustering(fsa_clustering, it->second, "results/" + directory + "/fsa_" + std::string(1, it->first) + ".txt");
        plot_clustering(cdba_clustering, it->second, "results/" + directory + "/cdba_" + std::string(1, it->first) + ".txt");
        plot_clustering(wedge_clustering, it->second, "results/" + directory + "/wedge_" + std::string(1, it->first) + ".txt");
        
        distance_t dba_cost = kMedianCost(curves, dba_clustering, integral_frechet<5>);
        distance_t fsa_cost = kMedianCost(curves, fsa_clustering, integral_frechet<5>);
        distance_t cdba_cost = kMedianCost(curves, cdba_clustering, integral_frechet<5>);
        distance_t wedge_cost = kMedianCost(curves, wedge_clustering, integral_frechet<5>);

        k_med_scores << it->first 
        << "," << std::to_string(dba_cost) 
        << "," << std::to_string(fsa_cost) 
        << "," << std::to_string(cdba_cost) 
        << "," << std::to_string(wedge_cost) 
        << "\n";

        plots << "dba_" + std::string(1, it->first) << std::endl;
        plots << "fsa_" + std::string(1, it->first) << std::endl;
        plots << "cdba_" + std::string(1, it->first) << std::endl;
        plots << "wedge_" + std::string(1, it->first) << std::endl;
    }
    k_med_scores.close();
    plots.close();
}

void curve_complexity_experiment_characters() {
    
    system("mkdir results");

    int n = 50;
    int k = 2;

    for (int l = 6; l <= 12; ++l) {
		std::string dir_name =
			"char_exp_" + std::to_string(l) + "_" + std::to_string(k) + "_" + std::to_string(n);
        center_update_experiment_characters(dir_name, n, k, l);
    }
}


void center_update_experiment_pigeons(std::string directory, int k, int l) {
    std::vector<std::string> pigeons = {
        "Bladon & Church route recapping/church hanborough/a94", "Bladon & Church route recapping/church hanborough/c22", "Bladon & Church route recapping/church hanborough/k77",
        "Bladon & Church route recapping/church hanborough/l29", "Bladon & Church route recapping/church hanborough/liv", "Bladon & Church route recapping/church hanborough/r47", "Bladon & Church route recapping/church hanborough/s93",
        "Horspath/H22", "Horspath/H27", "Horspath/H30", "Horspath/H35", "Horspath/H38", "Horspath/H41", "Horspath/H42", "Horspath/H71",
        "Weston/H23", "Weston/H31", "Weston/H32", "Weston/H34", "Weston/H36", "Weston/H50", "Weston/H58", "Weston/H62",
        "Bladon & Church route recapping/bladon heath/p39", "Bladon & Church route recapping/church hanborough/c70", "Bladon & Church route recapping/bladon heath/p29", "Bladon & Church route recapping/bladon heath/a55", "Bladon & Church route recapping/bladon heath/brc", "Bladon & Church route recapping/bladon heath/c17", "Bladon & Church route recapping/bladon heath/c35",
        "Bladon & Church route recapping/bladon heath/p94"
    };

    system(("mkdir results/" + directory).c_str());
     system(("mkdir 'results/" + directory + "/Bladon & Church route recapping'").c_str());
    system(("mkdir 'results/" + directory + "/Bladon & Church route recapping/church hanborough'").c_str());
    system(("mkdir 'results/" + directory + "/Bladon & Church route recapping/bladon heath'").c_str());
    system(("mkdir results/" + directory + "/Weston").c_str());
    system(("mkdir results/" + directory + "/Horspath").c_str());

    std::fstream clustering_names;
    clustering_names.open("results/" + directory + "/plots.txt", std::fstream::out | std::fstream::trunc);
    std::fstream k_med_scores;
    k_med_scores.open("results/" + directory + "/k_med_scores.txt", std::fstream::out | std::fstream::trunc);

    k_med_scores << "pigeons\tdba\tfsa\tcdba\twedge\n";

    for (std::string pigeon: pigeons) {
        Curves raw_curves = io::read_pigeon_curves_utm("data/Data_for_Mann_et_al_RSBL 2/" + pigeon);
        Curves curves = Curves();
        for (auto curve: raw_curves) {
            curves.push_back(curve.naive_l_simplification(10));
        }

        Clustering initial_clustering = computeClustering(curves, k, l, ClusterAlg::Pam, integral_frechet<250>, "", false);

        Clustering dba_clustering = initial_clustering;
        int count = 1;
        int max_count = 10;
        while (count <= max_count && computerCenters(curves, dba_clustering, l, CenterAlg::dba, dtw)) {
            updateClustering(curves, dba_clustering, dtw, nullptr);
            ++count;
        }

        Clustering fsa_clustering = initial_clustering;
        count = 1;
        while (count <= max_count && computerCenters(curves, fsa_clustering, l, CenterAlg::fCenter, frechet)) {
            updateClustering(curves, fsa_clustering, frechet, nullptr);
            ++count;
        }

        Clustering cdba_clustering = initial_clustering;
        count = 1;
        while (count <= max_count && computerCenters(curves, cdba_clustering, l, CenterAlg::cdbaPigeon, integral_frechet<250>)) {
            updateClustering(curves, cdba_clustering, integral_frechet<250>, nullptr);
            ++count;
        }

        Clustering wedge_clustering = initial_clustering;
        count = 1;
        while (count <= max_count && computerCenters(curves, wedge_clustering, l, CenterAlg::wedgePigeon, integral_frechet<250>)) {
            updateClustering(curves, wedge_clustering, integral_frechet<250>, nullptr);
            ++count;
        }


        std::cout << "results/" + directory + "/" +  pigeon + "_dba.txt" << "\n";
        plot_clustering(dba_clustering, curves, "results/" + directory + "/" +  pigeon + "_dba.txt");
        plot_clustering(fsa_clustering, curves, "results/" + directory + "/" +  pigeon + "_fsa.txt");
        plot_clustering(cdba_clustering, curves, "results/" + directory + "/" +  pigeon + "_cdba_integral.txt");
        plot_clustering(wedge_clustering, curves, "results/" + directory + "/" +  pigeon + "_wedge_integral.txt");

        clustering_names << pigeon + "_dba\n";
        clustering_names << pigeon + "_fsa\n";
        clustering_names << pigeon + "_cdba_integral\n";
        clustering_names << pigeon + "_wedge_integral\n";

        std::vector<Clustering> clusterings = {dba_clustering, fsa_clustering, cdba_clustering, wedge_clustering};

        k_med_scores << pigeon << "\t";
        for (int i = 0; i < clusterings.size(); ++i) {
            Clustering clustering = clusterings[i];
            distance_t k_med_score = 0;
            for (auto& cluster: clustering) {
                k_med_score += calcC2CDist(curves, cluster.center_curve, cluster.curve_ids, C2CDist::Median, integral_frechet<250>);
            }

            k_med_scores << std::to_string(k_med_score);
            
            if (i < clusterings.size() - 1) {
                k_med_scores << "\t";
            } else {
                k_med_scores << "\n";
            }
        }
    }

    k_med_scores.close();
    clustering_names.close();
}

void curve_complexity_experiment_pigeons() {
    
    system("mkdir results");

    int k = 3;
    
    for (int l = 6; l <= 12; ++l) {
        center_update_experiment_pigeons("pigeon_exp_" + std::to_string(l) + "_" + std::to_string(k), k, l);
    }

}


void running_time_experiment() {
    std::vector<std::string> pigeon_datasets = {
        "data/Data_for_Mann_et_al_RSBL 2/Bladon & Church route recapping/church hanborough/utm",
        "data/Data_for_Mann_et_al_RSBL 2/Bladon & Church route recapping/bladon heath/utm",
        "data/Data_for_Mann_et_al_RSBL 2/Horspath/utm",
        "data/Data_for_Mann_et_al_RSBL 2/Weston/utm"
    };

    std::vector<std::string> names = {
        "bladon_heath", "church_hanborough",
        "horspath", "weston"
    };


    std::fstream output_file;
    output_file.open("running_time.txt", std::fstream::out | std::fstream::trunc);

    for (int i = 0; i < pigeon_datasets.size(); ++i) {
        auto pigeon_data = pigeon_datasets[i];
        Curves raw_curves = io::read_pigeon_curves(pigeon_data);
        Curves curves = Curves();
        for (auto curve: raw_curves) {
            curves.push_back(curve.naive_l_simplification(100));
        }

        auto start_time = std::chrono::high_resolution_clock::now();
        Clustering clustering = computeCenterClustering(curves, 5, 10, ClusterAlg::Gonzalez, CenterAlg::cdba, integral_frechet<250>, integral_frechet<250>, "", 1);
        auto end_time = std::chrono::high_resolution_clock::now();
        auto total_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
        output_file << names[i] << "\t" << total_time << "\n";
    }

    output_file.close();
}
