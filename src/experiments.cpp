#include "experiments.h"
#include <fstream>
#include <map>
#include <stdexcept>
#include <vector>

#include "Curve.h"
#include "clustering/center_algs.h"
#include "clustering/center_clustering_algs.h"
#include "clustering/plot_clustering.h"
#include "distance_functions.h"
#include "simplification/agarwal.h"
#include "synthetic_curves.h"
#include "utils/io.h"
#include "utils/remove_stops.h"

namespace {
    distance_t kMedianCost(Curves const& curves, Clustering const& clustering,
            std::function<distance_t(Curve const&, Curve const&)> const& dist) {
        distance_t cost = 0.0;
        #pragma omp parallel for schedule(dynamic) reduction(+:cost)
        for (std::size_t i = 0; i < clustering.size(); ++i) {
            auto const& cluster = clustering[i];
            cost += clustering::calcC2CDist(curves, cluster.center_curve,
                cluster.curve_ids, clustering::C2CDist::Median, dist);
        }
        return cost;
    }

    distance_t kMedianCostMat(Clustering const& clustering,
        DistanceMatrix<distance_t> const& dist_matrix) {
        distance_t cost = 0.0;
        for (auto const& cluster: clustering)
            for (auto const& curve_id: cluster.curve_ids)
                cost += dist_matrix.at(curve_id, cluster.center_id);

        return cost;
    }

    void compute_simp_matrix(Curves const& curves, std::size_t l,
            std::function<distance_t(Curve const&, Curve const&)> const& dist,
            std::string const& out) {
        Curves simplifications;
        simplifications.reserve(curves.size());
        for (auto const& curve: curves) {
            simplifications.emplace_back(
                simplification::greedy::simplify(curve, l, df::frechet_lt));
        }
        CurveSimpMatrix matrix(curves, simplifications, dist);
        matrix.write(out);
    }

    [[maybe_unused]] void compute_dist_matrix(Curves const& curves,
            std::function<distance_t(Curve const&, Curve const&)> const& dist,
            std::string const& out) {
        SymmetricMatrix matrix(curves.size());
        #pragma omp parallel for schedule(dynamic)
        for (std::size_t i = 0; i < curves.size(); ++i) {
            matrix.at(i, i) = 0.0;
            for (std::size_t j = i + 1; j < curves.size(); ++j)
                matrix.at(i, j) = dist(curves[i], curves[j]);
        }
        matrix.write(out);
    }

    CurveSimpMatrix read_or_create(std::string const& file,
            Curves const& curves, std::size_t l,
            std::function<distance_t(Curve const&, Curve const&)> const& dist) {
        std::ifstream _file(file);
        if (!_file.is_open())
            compute_simp_matrix(curves, l, dist, file);
        else
            std::cout << "file found\n";
        return CurveSimpMatrix(file);
    }

    // SymmetricMatrix read_or_create_sym(std::string const& file,
    //         Curves const& curves, std::size_t l,
    //         std::function<distance_t(Curve const&, Curve const&)> const& dist) {
    //     // if (!std::filesystem::exists(file))
    //         compute_simp_matrix(curves, l, dist, file);
    //     return CurveSimpMatrix(file);
    // }

    [[maybe_unused]] SymmetricMatrix read_or_create(std::string const& file,
            Curves const& curves,
            std::function<distance_t(Curve const&, Curve const&)> const& dist) {
        compute_dist_matrix(curves, dist, file);
        std::ifstream matrix(file);
        if (!matrix.is_open())
            throw std::runtime_error("Failed to open file " + file);
        return SymmetricMatrix::read(matrix);
    }
}


distance_t compute_aic(Curves& curves, int k, int l, std::string matrix_path, 
bool fix_start, bool fix_end) {
    auto const dummy_dist = [](Curve const&, Curve const&) noexcept {
        return 0.0;
    };
    CurveSimpMatrix dist_matrix = read_or_create(matrix_path,
        curves, 10, df::integral_frechet);

    Clustering initial = clustering::computeClustering(curves, k, 10,
            clustering::ClusterAlg::PAM, dist_matrix, dummy_dist,
            df::frechet_lt);


    Clustering cdba_res = clustering::computeCenterClustering(curves, k, 10,
            fix_start, fix_end, clustering::ClusterAlg::PAM,
            clustering::CenterAlg::cdba, dist_matrix, dummy_dist,
            df::integral_frechet, df::frechet_lt, 1);


    return kMedianCost(curves, cdba_res, df::integral_frechet) + 40*k;
}

distance_t compute_ch_index(Curves& curves, int k, int l, std::string matrix_path,
bool fix_start, bool fix_end) {
    auto const dummy_dist = [](Curve const&, Curve const&) noexcept {
        return 0.0;
    };
    CurveSimpMatrix dist_matrix = read_or_create(matrix_path,
    curves, l, df::integral_frechet);

    Clustering initial = clustering::computeClustering(curves, k, 10,
            clustering::ClusterAlg::PAM, dist_matrix, dummy_dist,
            df::frechet_lt);

    Clustering cdba_res = clustering::computeCenterClustering(curves, k, 10,
            fix_start, fix_end, clustering::ClusterAlg::PAM,
            clustering::CenterAlg::cdba, dist_matrix, dummy_dist,
            df::integral_frechet, df::frechet_lt, 1);

    distance_t B = 0;
    distance_t W = 0;

    int B_count = 0;
    int W_count = 0;

    for (int i = 0; i < cdba_res.size(); ++i) {
        auto& cluster = cdba_res[i];
        for (int id_index = 0; id_index < cluster.curve_ids.size(); ++id_index) {
            
            for (int other_index = id_index + 1; other_index < cluster.curve_ids.size(); ++other_index) {
                W += df::integral_frechet(
                    curves[cluster.curve_ids[id_index]],
                    curves[cluster.curve_ids[other_index]]
                );
                W_count++;
            }

            CurveID id = cluster.curve_ids[id_index];

            for (int j = i+1; j < cdba_res.size(); ++j) {
                auto& other_cluster = cdba_res[j];
                for (CurveID o_id: other_cluster.curve_ids) {
                    Curve a = curves[id];
                    Curve b = curves[o_id];
                    B += df::integral_frechet(a, b);
                    B_count++;
                }
            }
        }
    }

    return (B / B_count) / (W / W_count);
    // return ((B / (k-1)) / (W / (curves.size()-k)));
}


void experiments::aic() {
    Curves base = io::read_curves("data/characters/data");

    Curves simplified = Curves();
    for (int i = 1500; i < 1560; ++i) {
        auto& curve = base[i];

        auto new_curve = Curve();
        for (auto p: curve.get_points())
            new_curve.push_back(p/15);

        simplified.push_back(new_curve.naive_l_simplification(50));
    }

    for (int i = 2000; i < 2060; ++i) {
        auto& curve = base[i];
        
        auto new_curve = Curve();
        for (auto p: curve.get_points())
            new_curve.push_back(p/15);

        simplified.push_back(new_curve.naive_l_simplification(50));
    }

    for (int i = 2500; i < 2560; ++i) {
        auto& curve = base[i];

        auto new_curve = Curve();
        for (auto p: curve.get_points())
            new_curve.push_back(p/15);

        simplified.push_back(new_curve.naive_l_simplification(50));
    }

    for (auto& c: simplified)
        std::cout << c.name() << std::endl;
    
    auto const dummy_dist = [](Curve const&, Curve const&) noexcept {
        return 0.0;
    };

    std::ofstream aic_results("elbow_graph_results/characters_aic.csv");
    aic_results << "k,aic\n";

    for (int k = 2; k <= 20; ++k) {
        std::cout << k << "\n";
        CurveSimpMatrix dist_matrix = read_or_create("matrices/3_60_characters_scaled.txt",
        simplified, 10, df::integral_frechet);

        Clustering initial = clustering::computeClustering(simplified, k, 10,
                clustering::ClusterAlg::PAM, dist_matrix, dummy_dist,
                df::frechet_lt);

        Clustering cdba_res = clustering::computeCenterClustering(simplified, k, 10,
                true, false, clustering::ClusterAlg::PAM,
                clustering::CenterAlg::cdba, dist_matrix, dummy_dist,
                df::integral_frechet, df::frechet_lt, 1);

        aic_results << k << "," << kMedianCost(simplified, cdba_res, df::integral_frechet) + 40*k << std::endl;
    }

    aic_results.close();
}


void experiments::elbow_graph() {
    Curves base = io::read_curves("data/characters/data");

    Curves simplified = Curves();
    for (int i = 1500; i < 1530; ++i) {
        auto& curve = base[i];
        simplified.push_back(curve.naive_l_simplification(50));
    }

    for (int i = 2000; i < 2030; ++i) {
        auto& curve = base[i];
        simplified.push_back(curve.naive_l_simplification(50));
    }

    for (int i = 2500; i < 2530; ++i) {
        auto& curve = base[i];
        simplified.push_back(curve.naive_l_simplification(50));
    }

    for (auto& c: simplified)
        std::cout << c.name() << std::endl;
    
    auto const dummy_dist = [](Curve const&, Curve const&) noexcept {
        return 0.0;
    };

    

    // int k = 1;


    // std::ofstream results("elbow_graph_results/a_characters_l.csv");

    // results << "l,k_med_score\n";


    // for (int ell = 5; ell <= 20; ++ell) {

    //     CurveSimpMatrix dist_matrix = read_or_create("matrices/a_characters_"+std::to_string(ell)+".txt",
    //     simplified, ell, df::integral_frechet);

    //     Clustering initial = clustering::computeClustering(simplified, k, ell,
    //             clustering::ClusterAlg::PAM, dist_matrix, dummy_dist,
    //             df::frechet_lt);

    //     Clustering cdba_res = clustering::computeCenterClustering(simplified, k, ell,
    //             false, false, clustering::ClusterAlg::PAM,
    //             clustering::CenterAlg::cdba, dist_matrix, dummy_dist,
    //             df::integral_frechet, df::frechet_lt, 1);

    //     results << ell << "," << kMedianCost(simplified, cdba_res, df::integral_frechet) << "\n";
    // }

    // results.close();


    // std::ofstream k_results("elbow_graph_results/a_characters_k.csv");

    // k_results << "k,k_med_score\n";

    // for (int k = 1; k <= 20; ++k) {
    //     CurveSimpMatrix dist_matrix = read_or_create("matrices/a_characters_10.txt",
    //     simplified, 10, df::integral_frechet);

    //     Clustering initial = clustering::computeClustering(simplified, k, 10,
    //             clustering::ClusterAlg::PAM, dist_matrix, dummy_dist,
    //             df::frechet_lt);

    //     Clustering cdba_res = clustering::computeCenterClustering(simplified, k, 10,
    //             false, false, clustering::ClusterAlg::PAM,
    //             clustering::CenterAlg::cdba, dist_matrix, dummy_dist,
    //             df::integral_frechet, df::frechet_lt, 1);

    //     k_results << k << "," << kMedianCost(simplified, cdba_res, df::integral_frechet) << "\n";
    // }

    // k_results.close();

    std::ofstream ch_results("elbow_graph_results/characters_ch.csv");


    ch_results << "k,ch_index\n";

    for (int k = 2; k <= 20; ++k) {
        std::cout << k << "\n";
        CurveSimpMatrix dist_matrix = read_or_create("matrices/3_30_characters.txt",
        simplified, 10, df::integral_frechet);

        Clustering initial = clustering::computeClustering(simplified, k, 10,
                clustering::ClusterAlg::PAM, dist_matrix, dummy_dist,
                df::frechet_lt);

        Clustering cdba_res = clustering::computeCenterClustering(simplified, k, 10,
                true, false, clustering::ClusterAlg::PAM,
                clustering::CenterAlg::cdba, dist_matrix, dummy_dist,
                df::integral_frechet, df::frechet_lt, 1);

        distance_t ch_index = 0;

        distance_t B = 0;
        distance_t W = 0;

        for (int i = 0; i < cdba_res.size(); ++i) {
            auto& cluster = cdba_res[i];
            for (int id_index = 0; id_index < cluster.curve_ids.size(); ++id_index) {
                
                for (int other_index = id_index + 1; other_index < cluster.curve_ids.size(); ++other_index) {
                    W += df::integral_frechet(
                        simplified[cluster.curve_ids[id_index]],
                        simplified[cluster.curve_ids[other_index]]
                    );
                }

                CurveID id = cluster.curve_ids[id_index];

                for (int j = i+1; j < cdba_res.size(); ++j) {
                    auto& other_cluster = cdba_res[j];
                    for (CurveID o_id: other_cluster.curve_ids) {
                        Curve a = simplified[id];
                        Curve b = simplified[o_id];
                        B += df::integral_frechet(a, b);
                    }
                }
            }
        }

        clustering::plot_clustering(cdba_res, simplified, "clustering_"+std::to_string(k)+".txt");

        ch_results << k << "," << ((B / (k-1)) / (W / (simplified.size()-k))) << "\n";
    }

    ch_results.close();
}


std::vector<std::vector<char>> 
find_subsets(std::vector<char> characters, int set_size, int set_count, int seed) {

    std::uniform_real_distribution<double> dis(0.0, 1.0);
    std::mt19937 generator (seed);

    std::vector<std::vector<char>> output = std::vector<std::vector<char>>();
    std::vector<int> indicator(characters.size(), 0);
    for (int i = 0; i < set_size; ++i)
        indicator[i] = 1;
    for (int i = 0; i < set_count; ++i) {
        std::vector<char> subset = std::vector<char>();
        std::shuffle(indicator.begin(), indicator.end(), generator);
        for (int j = 0; j < indicator.size(); ++j) {
            int ind = indicator[j];
            if (ind == 1) {
                subset.push_back(characters[j]);
            }
        }
        output.push_back(subset);
    }

    return output;
}


void experiments::characters_param_selection() {
    int n = 100;
    Curves curves = io::read_curves("data/characters/data");
    std::map<char, Curves> curves_by_letter;
    for (auto const& curve: curves) {
        char const letter = curve.name()[21];
        if (curves_by_letter[letter].size() < n)
            curves_by_letter[letter].emplace_back(curve.naive_l_simplification(30));
    }

    std::vector<char> characters = std::vector<char>();

    for (auto& pair: curves_by_letter) {
        characters.push_back(pair.first);
    }

    int seed = 10;

    for (int set_size = 3; set_size <= 3; set_size++) {
        // int set_count = 11 - set_size;
        int set_count = 4;
        auto subsets = find_subsets(characters, set_size, set_count, seed);

        for (auto& set: subsets) {
            Curves curve_subset = Curves();
            std::string set_string = "";
            for (auto ch: set) {
                std::cout << ch << std::endl;
                std::string ch_s(1, ch);
                set_string += ch_s;
                set_string += "_";
                for (int i = 0; i < 20; ++i) {
                    curve_subset.push_back(curves_by_letter[ch][i]);
                }
            }

            std::cout << set_string << std::endl;

            std::ofstream results("param_selection/characters_" + set_string + ".txt");
            results << "k,aic,ch_index\n";
            for (int k = 3; k <= 20; ++k) {
                std::cout << k << std::endl;
                std::cout << curve_subset.size() << std::endl;
                std::string matrix_path = "matrices/" + set_string + ".txt";
                auto ch_index = compute_ch_index(curve_subset, k, 10, matrix_path, true, false);
                std::cout << "computed ch_index...\n";
                auto aic = compute_aic(curve_subset, k, 10, matrix_path, true, false);
                std::cout << "computed aic...\n";
                results << k << "," << aic << "," << ch_index << std::endl;
            }
            results.close();
        }
    }
}