#include "experiments.h"

#include <filesystem>
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
        if (!std::filesystem::exists(file))
            compute_simp_matrix(curves, l, dist, file);
        return CurveSimpMatrix(file);
    }

    [[maybe_unused]] SymmetricMatrix read_or_create(std::string const& file,
            Curves const& curves,
            std::function<distance_t(Curve const&, Curve const&)> const& dist) {
        if (!std::filesystem::exists(file))
            compute_dist_matrix(curves, dist, file);
        std::ifstream matrix(file);
        if (!matrix.is_open())
            throw std::runtime_error("Failed to open file " + file);
        return SymmetricMatrix::read(matrix);
    }
}

void experiments::synthetic_curve_experiment() {
    std::filesystem::create_directory("synthetic_curves");
    if (!std::filesystem::exists("synthetic_curves/dataset.txt")) {
        Curve base = io::read_curves("data/characters/data")[0];
        synth::write_curves(base.naive_l_simplification(200));
    }
    Curves curves = io::read_curves("synthetic_curves");
    CurveSimpMatrix dummy_matrix({}, {},
        [](Curve const&, Curve const&) noexcept { return 0.0; });
    Clustering gonzalez_clustering = clustering::computeCenterClustering(
        curves, 1, 15, false, false, clustering::ClusterAlg::Gonzalez,
        clustering::CenterAlg::cdba, dummy_matrix,
        df::integral_frechet_fast, df::integral_frechet, df::frechet_lt, 1);
    clustering::plot_clustering(gonzalez_clustering, curves, "plot.txt");
}

void experiments::initial_clustering_experiment() {
    std::array pigeon_sites {
        "Bladon & Church route recapping/bladon heath",
        "Bladon & Church route recapping/church hanborough",
        "Horspath",
        "Weston"
    };

    std::array dist_matrix_paths {
        "bladon_heath_10.txt",
        "church_10.txt",
        "horspath_10.txt",
        "weston_10.txt"
    };

    std::array l1_matrix_paths {
        "bladon_heath_l1_10.txt",
        "church_l1_10.txt",
        "horspath_l1_10.txt",
        "weston_l1_10.txt"
    };

    std::filesystem::create_directories(
        "results/initial_clustering_experiment/figures");

    std::string const filename = 
        "results/initial_clustering_experiment/results.txt";
    std::ofstream results(filename);

    std::string const l1_filename = 
        "results/initial_clustering_experiment/results_l1.txt";
    std::ofstream l1_results(l1_filename);
    if (!results.is_open())
        throw std::runtime_error("Failed to open file " + filename);

    if (!l1_results.is_open())
        throw std::runtime_error("Failed to open file " + filename);

    for (std::size_t i = 0 ; i < pigeon_sites.size(); ++i) {
        std::string const& site = pigeon_sites[i];
        std::string const& dist_matrix_path = dist_matrix_paths[i];
        std::string const& l1_matrix_path = l1_matrix_paths[i];

        Curves raw_pigeon_curves = io::read_curves(
            "data/Data_for_Mann_et_al_RSBL/" + site + "/utm", 1);
        Curves curves;
        curves.reserve(raw_pigeon_curves.size());
        for (auto const& curve: raw_pigeon_curves)
            curves.emplace_back(curve.naive_l_simplification(200));

        auto const dummy_dist = [](Curve const&, Curve const&) noexcept {
            return 0.0;
        };
        CurveSimpMatrix dist_matrix = read_or_create(dist_matrix_path, curves,
            10, df::integral_frechet);

        CurveSimpMatrix l1_matrix = read_or_create(l1_matrix_path, curves,
        10, df::heur_cdtw_2d_l1_l1);

        for (std::size_t k = 1; k <= 10; ++k) {
            // Run the experiment 6 times without center updates, helps with
            // randomness in Gonzalez initialization.
            Clustering gonzalez_res =
                clustering::computeCenterClustering(curves, k, 10, true, true,
                    clustering::ClusterAlg::Gonzalez,
                    clustering::CenterAlg::none, dist_matrix,
                    df::integral_frechet, df::integral_frechet,
                    df::frechet_lt, 6);

            Clustering gonzalez_pam_res =
                clustering::computeCenterClustering(curves, k, 10, true, true,
                    clustering::ClusterAlg::GonzalezPAM,
                    clustering::CenterAlg::none, dist_matrix,
                    df::integral_frechet, df::integral_frechet,
                    df::frechet_lt, 6);
            // No randomness in PAM, so running once is enough.
            Clustering pam_res =
                clustering::computeClustering(curves, k, 10,
                    clustering::ClusterAlg::PAM, dist_matrix, dummy_dist,
                    df::frechet_lt);

            // Run the experiment 6 times without center updates, helps with
            // randomness in Gonzalez initialization.
            Clustering gonzalez_res_l1 =
                clustering::computeCenterClustering(curves, k, 10, true, true,
                    clustering::ClusterAlg::Gonzalez,
                    clustering::CenterAlg::none, l1_matrix,
                    df::heur_cdtw_2d_l1_l1, df::heur_cdtw_2d_l1_l1,
                    df::frechet_lt, 6);

            Clustering gonzalez_pam_res_l1 =
                clustering::computeCenterClustering(curves, k, 10, true, true,
                    clustering::ClusterAlg::GonzalezPAM,
                    clustering::CenterAlg::none, l1_matrix,
                    df::heur_cdtw_2d_l1_l1, df::heur_cdtw_2d_l1_l1,
                    df::frechet_lt, 6);
            // No randomness in PAM, so running once is enough.
            Clustering pam_res_l1 =
                clustering::computeClustering(curves, k, 10,
                    clustering::ClusterAlg::PAM, l1_matrix, dummy_dist,
                    df::frechet_lt);


            results << site << "\t" << k << "\t" <<
            kMedianCostMat(gonzalez_pam_res, dist_matrix) <<
            "\t" << kMedianCostMat(pam_res, dist_matrix) <<
            "\t" << kMedianCostMat(gonzalez_res, dist_matrix) << std::endl;

            l1_results << site << "\t" << k << "\t" <<
            kMedianCostMat(gonzalez_pam_res_l1, l1_matrix) <<
            "\t" << kMedianCostMat(pam_res_l1, l1_matrix) <<
            "\t" << kMedianCostMat(gonzalez_res_l1, l1_matrix) << std::endl;
        }
    }
}

void experiments::center_update_experiment_movebank(
        std::string const& directory, std::size_t k, std::size_t l,
        bool remove_stops) {
    Curves raw_curves = io::read_curves("data/movebank/projection");
    Curves curves;
    curves.reserve(raw_curves.size());

    // Simplify and scale trajectories.
    for (auto const& curve: raw_curves) {
        auto simplified = curve.naive_l_simplification(400);
        Points transformed;
        transformed.reserve(curve.size());
        for (auto const& p: simplified.get_points())
            transformed.emplace_back(p.x / 1000, p.y / 1000);
        if (remove_stops)
            curves.emplace_back(io::remove_stops(Curve(transformed)));
        else
            curves.emplace_back(transformed);
    }

    std::string const path = "results/" + directory + "/";
    std::filesystem::create_directories(path + "figures");
    std::string const filename(path + "plots.txt");
    std::ofstream clustering_names(filename);
    if (!clustering_names.is_open())
        throw std::runtime_error("Failed to open file " + filename);
    clustering_names << "cdba\nwedge\ndba\nfsa";

    auto const dummy_dist = [](Curve const&, Curve const&) noexcept {
        return 0.0;
    };
    CurveSimpMatrix dist_matrix = read_or_create("movebank.txt", curves, l,
        df::integral_frechet);

    CurveSimpMatrix l1_matrix = read_or_create("movebank_l1.txt", curves, l,
        df::heur_cdtw_2d_l1_l1);

    Clustering initial = clustering::computeClustering(curves, k, l,
        clustering::ClusterAlg::PAM, dist_matrix, dummy_dist, df::frechet_lt);

    // PAM is not random, so we can just recompute the initial clustering
    // four times.
    Clustering cdba_res = clustering::computeCenterClustering(curves, k, l,
        false, true, clustering::ClusterAlg::PAM, clustering::CenterAlg::cdba,
        dist_matrix, dummy_dist, df::integral_frechet, df::frechet_lt, 1);

    Clustering cdba_l1_res = clustering::computeCenterClustering(curves, k, l,
        false, true, clustering::ClusterAlg::PAM, clustering::CenterAlg::cdba_l1,
        l1_matrix, dummy_dist, df::heur_cdtw_2d_l1_l1, df::frechet_lt, 1);

    Clustering dba_res = clustering::computeCenterClustering(curves, k, l,
        false, true, clustering::ClusterAlg::PAM, clustering::CenterAlg::dba,
        dist_matrix, dummy_dist, df::dtw, df::frechet_lt, 1);

    Clustering fsa_res = clustering::computeCenterClustering(curves, k, l,
        false, true, clustering::ClusterAlg::PAM, clustering::CenterAlg::fsa,
        dist_matrix, dummy_dist, df::frechet, df::frechet_lt, 1);

    Clustering wedge_res = clustering::computeCenterClustering(curves, k, l,
        false, true, clustering::ClusterAlg::PAM, clustering::CenterAlg::wedge,
        dist_matrix, dummy_dist, df::integral_frechet, df::frechet_lt, 1,
        5, 20);

    clustering::plot_clustering(initial, curves, path + "initial.txt");
    clustering::plot_clustering(cdba_res, curves, path + "cdba.txt");
    clustering::plot_clustering(cdba_l1_res, curves, path + "cdba_l1.txt");
    clustering::plot_clustering(dba_res, curves, path + "dba.txt");
    clustering::plot_clustering(fsa_res, curves, path + "fsa.txt");
    clustering::plot_clustering(wedge_res, curves, path + "wedge.txt");
}


void experiments::center_update_experiment_characters(
        std::string const& directory, std::size_t n, std::size_t k,
        std::size_t l) {
    Curves curves = io::read_curves("data/characters/data");
    std::map<char, Curves> curves_by_letter;
    for (auto const& curve: curves) {
        char const letter = curve.name()[21];
        if (curves_by_letter[letter].size() < n)
            curves_by_letter[letter].emplace_back(curve);
    }

    std::string const path = "results/" + directory + "/";
    std::filesystem::create_directories(path + "figures");
    std::ofstream plots(path + "plots.txt");
    if (!plots.is_open())
        throw std::runtime_error("Failed to open file " + path + "plots.txt");

    std::ofstream k_med_scores(path + "k_med_scores.csv");
    if (!k_med_scores.is_open())
        throw std::runtime_error("Failed to open file " + path +
            "k_med_scores.csv");
    k_med_scores << "characters,DBA,FSA,CDBA,Wedge,initial\n";

    auto const dummy_dist = [](Curve const&, Curve const&) noexcept {
        return 0.0;
    };
    std::filesystem::create_directory("characters");
    for (auto const& [letter, crv]: curves_by_letter) {
        CurveSimpMatrix dist_matrix = read_or_create(std::string("characters/")
            + letter + ".txt", crv, l, df::integral_frechet);

        CurveSimpMatrix l1_matrix = read_or_create(std::string("characters/")
            + letter + ".l1_txt", crv, l, df::heur_cdtw_2d_l1_l1);
        Clustering initial = clustering::computeClustering(crv, k, l,
            clustering::ClusterAlg::PAM, dist_matrix, dummy_dist,
            df::frechet_lt);

        // PAM is not random, so we can just recompute the initial clustering
        // four times.
        Clustering cdba_res = clustering::computeCenterClustering(crv, k, l,
            false, false, clustering::ClusterAlg::PAM,
            clustering::CenterAlg::cdba, dist_matrix, dummy_dist,
            df::integral_frechet, df::frechet_lt, 1);
        
        Clustering cdba_l1_res = clustering::computeCenterClustering(crv, k, l,
            false, false, clustering::ClusterAlg::PAM,
            clustering::CenterAlg::cdba_l1, l1_matrix, dummy_dist,
            df::heur_cdtw_2d_l1_l1, df::frechet_lt, 1);


        Clustering dba_res = clustering::computeCenterClustering(crv, k, l,
            false, false, clustering::ClusterAlg::PAM,
            clustering::CenterAlg::dba, dist_matrix, dummy_dist,
            df::dtw, df::frechet_lt, 1);

        Clustering fsa_res = clustering::computeCenterClustering(crv, k, l,
            false, false, clustering::ClusterAlg::PAM,
            clustering::CenterAlg::fsa, dist_matrix, dummy_dist,
            df::frechet, df::frechet_lt, 1);

        Clustering wedge_res = clustering::computeCenterClustering(crv, k, l,
            false, false, clustering::ClusterAlg::PAM,
            clustering::CenterAlg::wedge, dist_matrix, dummy_dist,
            df::integral_frechet, df::frechet_lt, 1, 0.125, 20);

        clustering::plot_clustering(dba_res, crv,
            path + "dba_" + letter + ".txt");
        clustering::plot_clustering(fsa_res, crv,
            path + "fsa_" + letter + ".txt");
        clustering::plot_clustering(cdba_res, crv,
            path + "cdba_" + letter + ".txt");
        clustering::plot_clustering(cdba_l1_res, crv,
            path + "cdba_l1_" + letter + ".txt");
        clustering::plot_clustering(wedge_res, crv,
            path + "wedge_" + letter + ".txt");
        clustering::plot_clustering(initial, crv,
            path + "init_" + letter + ".txt");

        k_med_scores << letter << ","
            << kMedianCost(crv, dba_res, df::integral_frechet) << ","
            << kMedianCost(crv, fsa_res, df::integral_frechet) << ","
            << kMedianCost(crv, cdba_res, df::integral_frechet) << ","
            << kMedianCost(crv, cdba_l1_res, df::integral_frechet) << ","
            << kMedianCost(crv, wedge_res, df::integral_frechet) << ","
            << kMedianCostMat(initial, dist_matrix) << std::endl;

        plots << "dba_" << letter << "\n" << "fsa_" << letter << "\n"
            << "cdba_" << letter << "\n" << "cdba_l1_" << letter << "\n"  << "wedge_" << letter << "\n"
            << "init_" << letter << std::endl;
    }
}

void experiments::curve_complexity_experiment_characters() {
    std::filesystem::create_directory("results");
    std::size_t n = 50, k = 2;
    for (std::size_t l = 6; l <= 12; ++l) {
        std::string const dir_name = "char_exp_" + std::to_string(l) + "_" +
            std::to_string(k) + "_" + std::to_string(n);
        center_update_experiment_characters(dir_name, n, k, l);
    }
}

void experiments::center_update_experiment_pigeons(
        std::string const& directory, std::size_t k, std::size_t l) {
    std::vector<std::string> p_bl {
        "a55", "brc", "c17", "c35", "p29", "p39", "p94"};
    std::vector<std::string> p_ch {
        "a94", "c22", "c70", "k77", "l29", "liv", "r47", "s93"};
    std::vector<std::string> p_hp {
        "H22", "H27", "H30", "H35", "H38", "H41", "H42", "H71"};
    std::vector<std::string> p_ws {
        "H23", "H31", "H32", "H34", "H36", "H50", "H58", "H62"};
    std::vector<std::vector<std::string>> sites {
        p_bl, p_ch, p_hp, p_ws};
    std::array<std::string, 4> site_paths {
        "Bladon & Church route recapping/bladon heath",
        "Bladon & Church route recapping/church hanborough",
        "Horspath", "Weston"};
    std::vector<std::string> pigeons;
    pigeons.reserve(p_bl.size() + p_ch.size() + p_hp.size() + p_ws.size());
    for (std::size_t i = 0; i < sites.size(); ++i) {
        for (std::size_t j = 0; j < sites[i].size(); ++j)
            pigeons.emplace_back(site_paths[i] + "/" + sites[i][j]);
        std::filesystem::create_directories(
            "results/" + directory + "/" + site_paths[i]);
    }

    std::ofstream plots("results/" + directory + "/plots.txt");
    if (!plots.is_open())
        throw std::runtime_error("Failed to open file results/" + directory +
            "plots.txt");
    std::ofstream k_med_scores("results/" + directory + "/k_med_scores.txt");
    if (!k_med_scores.is_open())
        throw std::runtime_error("Failed to open file results/" + directory +
            "k_med_scores.csv");
    k_med_scores << "pigeons\tdba\tfsa\tcdba\twedge\tinitial\n";

    auto const dummy_dist = [](Curve const&, Curve const&) noexcept {
        return 0.0;
    };
    for (auto const& site: site_paths)
        std::filesystem::create_directories("pigeons/" + site);
    for (std::string const& pigeon: pigeons) {
        Curves raw_curves = io::read_curves("data/Data_for_Mann_et_al_RSBL/" +
            pigeon + "/utm", 1);
        Curves curves;
        curves.reserve(raw_curves.size());
        for (auto const& curve: raw_curves)
            curves.emplace_back(curve.naive_l_simplification(200));

        CurveSimpMatrix dist_matrix = read_or_create(
            "pigeons/" + pigeon + ".txt", curves, l, df::integral_frechet);
        Clustering initial = clustering::computeClustering(curves, k, l,
            clustering::ClusterAlg::PAM, dist_matrix, dummy_dist,
            df::frechet_lt);


        CurveSimpMatrix l1_matrix = read_or_create(
            "pigeons/" + pigeon + ".l1_txt", curves, l, df::heur_cdtw_2d_l1_l1);

        // PAM is not random, so we can just recompute the initial clustering
        // four times.
        Clustering cdba_res = clustering::computeCenterClustering(curves, k, l,
            true, true, clustering::ClusterAlg::PAM,
            clustering::CenterAlg::cdba, dist_matrix, dummy_dist,
            df::integral_frechet, df::frechet_lt, 1);

        Clustering cdba_l1_res = clustering::computeCenterClustering(curves, k, l,
            true, true, clustering::ClusterAlg::PAM,
            clustering::CenterAlg::cdba_l1, l1_matrix, dummy_dist,
            df::heur_cdtw_2d_l1_l1, df::frechet_lt, 1);

        Clustering dba_res = clustering::computeCenterClustering(curves, k, l,
            true, true, clustering::ClusterAlg::PAM,
            clustering::CenterAlg::dba, dist_matrix, dummy_dist,
            df::dtw, df::frechet_lt, 1);

        Clustering fsa_res = clustering::computeCenterClustering(curves, k, l,
            true, true, clustering::ClusterAlg::PAM,
            clustering::CenterAlg::fsa, dist_matrix, dummy_dist,
            df::frechet, df::frechet_lt, 1);

        Clustering wedge_res = clustering::computeCenterClustering(curves, k, l,
            true, true, clustering::ClusterAlg::PAM,
            clustering::CenterAlg::wedge, dist_matrix, dummy_dist,
            df::integral_frechet, df::frechet_lt, 1, 5, 10);

        clustering::plot_clustering(dba_res, curves,
            "results/" + directory + "/" + pigeon + "_dba.txt");
        clustering::plot_clustering(fsa_res, curves,
            "results/" + directory + "/" + pigeon + "_fsa.txt");
        clustering::plot_clustering(cdba_res, curves,
            "results/" + directory + "/" + pigeon + "_cdba.txt");
        clustering::plot_clustering(cdba_l1_res, curves,
            "results/" + directory + "/" + pigeon + "_cdba_l1.txt");
        clustering::plot_clustering(wedge_res, curves,
            "results/" + directory + "/" + pigeon + "_wedge.txt");
        clustering::plot_clustering(initial, curves,
            "results/" + directory + "/" + pigeon + "_init.txt");

        plots << pigeon << "_dba\n" << pigeon << "_fsa\n"
            << pigeon << "_cdba\n" << pigeon << "_cdba_l1\n" << pigeon << "_wedge\n"
            << pigeon <<"_init" << std::endl;

        k_med_scores << pigeon << "\t"
            << kMedianCost(curves, dba_res, df::integral_frechet) << "\t"
            << kMedianCost(curves, fsa_res, df::integral_frechet) << "\t"
            << kMedianCost(curves, cdba_res, df::integral_frechet) << "\t"
            << kMedianCost(curves, cdba_l1_res, df::integral_frechet) << "\t"
            << kMedianCost(curves, wedge_res, df::integral_frechet) << "\t"
            << kMedianCostMat(initial, dist_matrix) << std::endl;
    }
}

void experiments::curve_complexity_experiment_pigeons() {
    std::size_t k = 3;
    for (std::size_t l = 6; l <= 12; l += 6)
        center_update_experiment_pigeons("pigeon_exp_" + std::to_string(l) +
            "_" + std::to_string(k), k, l);
}

void experiments::find_wedge_params_pigeons() {
    std::string const pigeon =
        "Bladon & Church route recapping/bladon heath/a55";
    Curves raw_curves = io::read_curves("data/Data_for_Mann_et_al_RSBL/" +
        pigeon + "/utm", 1);
    Curves curves;
    curves.reserve(raw_curves.size());
    for (auto const& curve: raw_curves)
        curves.emplace_back(curve.naive_l_simplification(200));

    std::size_t k = 3, l = 10;
    auto const dummy_dist = [](Curve const&, Curve const&) noexcept {
        return 0.0;
    };
    std::filesystem::create_directory("pigeons");
    CurveSimpMatrix dist_matrix = read_or_create("pigeons/" + pigeon + ".txt",
        curves, l, df::integral_frechet);
    Clustering initial_clustering = clustering::computeClustering(curves, k, l,
        clustering::ClusterAlg::PAM, dist_matrix, dummy_dist, df::frechet_lt);
    auto baseline = kMedianCostMat(initial_clustering, dist_matrix);
    distance_t best_eps = 0.0;
    int best_r = 0;

    std::array<distance_t, 4> epsilons {10, 30, 50, 70};
    std::array<int, 5> radii {5, 10, 20, 30, 40};
    #pragma omp parallel for collapse(2) schedule(dynamic)
    for (unsigned i = 0; i < epsilons.size(); ++i) {
        for (unsigned j = 0; j < radii.size(); ++j) {
            distance_t eps = epsilons[i];
            int radius = radii[j];
            Clustering start(initial_clustering);
            bool res = clustering::computeCenters(curves, start,
                clustering::CenterAlg::wedge, true, true, df::integral_frechet,
                eps, radius);
            if (!res)
                continue;
            auto new_cost = kMedianCost(curves, start, df::integral_frechet);
            #pragma omp critical(find_best_wedge)
            if (new_cost < baseline) {
                baseline = new_cost;
                best_eps = eps;
                best_r = radius;
            }
        }
    }
    std::cout << "Best improvement for pigeon a55 with eps = " +
        std::to_string(best_eps) + ", r = " + std::to_string(best_r)  << "."
        << std::endl;
}
