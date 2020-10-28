#include "classification_experiment.h"

#include <algorithm>
#include <chrono>
#include <iostream>
#include <random>
#include <unordered_map>
#include <vector>

#include "clustering/center_clustering_algs.h"
#include "utils/io.h"
#include "distance_functions.h"

#include "utils/defs.h"

using hrc = std::chrono::high_resolution_clock;
using time_point = hrc::time_point;
using ns = std::chrono::nanoseconds;
using ms = std::chrono::milliseconds;

Curves classification::sample(Curves const& curves, unsigned period) {
    Curves samples;
    for (std::size_t i = 0; i < curves.size(); i += period)
        samples.push_back(curves[i]);
    return samples;
}

void classification::classify_characters(
        clustering::ClusterAlg cluster_alg, clustering::CenterAlg center_alg,
        std::function<distance_t(Curve const&, Curve const&)> const& init_dist,
        std::function<distance_t(Curve const&, Curve const&)> const& dist) {
    unsigned const k = 6;
    unsigned const l = 12;
    unsigned const k_xv = 5; // the k of the cross-validation

    auto curves = sample(io::read_curves("data/characters/data"), 15);
    DEBUG("Loaded curves");

    std::random_device rd;
    std::mt19937_64 g(rd());
    std::shuffle(curves.begin(), curves.end(), g);

    // Put curves into own vectors by character.
    std::unordered_map<char, Curves> char_to_curves;
    for (auto const& curve: curves)
        char_to_curves[curve.name().front()].push_back(curve);

    // Upper level: characters, second level: CV folds
    std::vector<std::vector<Curves>> curves_by_char;
    for (auto const& pair: char_to_curves) {
        curves_by_char.emplace_back(k_xv);
        auto& cvs = pair.second;
        for (CurveID i = 0; i < cvs.size(); ++i)
            curves_by_char.back()[i%k_xv].push_back(std::move(cvs[i]));
    }

    unsigned correct = 0;
    unsigned wrong = 0;
    std::chrono::nanoseconds time(0);
    for (unsigned i = 0; i < k_xv; ++i) {
        Curves centers;
        for (auto const& curves_slots: curves_by_char) {
            Curves training_curves;
            for (unsigned j = 0; j < k_xv; ++j) {
                if (i == j)
                    continue;
                training_curves.insert(training_curves.end(),
                    curves_slots[j].begin(), curves_slots[j].end());
            }
            auto clustering = clustering::computeCenterClustering(
                training_curves, k, l, false, false, cluster_alg, center_alg,
                CurveSimpMatrix({}, {}, df::frechet), init_dist, dist,
                df::frechet_lt, 1);
            for (auto& cluster: clustering)
                centers.push_back(std::move(cluster.center_curve));

            // for simplicity, fill up centers using the first such that there
            // are k in total.
            if (clustering.empty())
                ERROR("The clustering has to contain at least one center.");

            auto representative = centers.back();
            while (centers.size() % k != 0)
                centers.push_back(representative);
        }

        // cross-validate on remaining data
        auto start = hrc::now();
        for (std::size_t char_id = 0; char_id < curves_by_char.size(); ++char_id) {
            auto const& cvs = curves_by_char[char_id][i];
            for (auto const& curve: cvs) {
                // classify curve
                distance_t min_dist = std::numeric_limits<distance_t>::max();
                std::size_t min_id = 0;
                for (CurveID center_id = 0; center_id < centers.size(); ++center_id) {
                    auto const& center = centers[center_id];
                    distance_t d = dist(curve, center);
                    if (d < min_dist) {
                        min_dist = d;
                        min_id = center_id;
                    }
                }

                // check if classification is correct
                if (char_id == min_id / k)
                    ++correct;
                else
                    ++wrong;
            }
        }
        time += std::chrono::duration_cast<ns>(hrc::now() - start);
    }

    std::cout << "Time: " << std::chrono::duration_cast<ms>(time).count()
              << " ms\nCorrect classifications: " << correct << "\n"
              << "Wrong classifications: " << wrong << "\nAccuracy: "
              << static_cast<double>(correct) / (correct + wrong) << std::endl;
}

void classification::identify_pigeons(
        clustering::ClusterAlg cluster_alg, clustering::CenterAlg center_alg,
        std::function<distance_t(Curve const&, Curve const&)> const& init_dist,
        std::function<distance_t(Curve const&, Curve const&)> const& dist) {

    std::vector<std::string> const p_bl {
        "a55", "brc", "c17", "c35", "p29", "p39", "p94"};
    std::vector<std::string> const p_ch {
        "a94", "c22", "c70", "k77", "l29", "liv", "r47", "s93"};
    std::vector<std::string> const p_hp {
        "H22", "H27", "H30", "H35", "H38", "H41", "H42", "H71"};
    std::vector<std::string> const p_ws {
        "H23", "H31", "H32", "H34", "H36", "H50", "H58", "H62"};
    std::vector<std::vector<std::string>> const sites {
        p_bl, p_ch, p_hp, p_ws};
    std::array<std::string, 4> const site_paths {
        "Bladon & Church route recapping/bladon heath",
        "Bladon & Church route recapping/church hanborough",
        "Horspath", "Weston"};

    std::vector<std::size_t> const k_bl {4, 3, 3, 4, 4, 5, 4};
    std::vector<std::size_t> const k_ch {4, 3, 5, 4, 4, 3, 4, 5};
    std::vector<std::size_t> const k_hp {5, 3, 4, 4, 4, 6, 4, 5};
    std::vector<std::size_t> const k_ws {6, 3, 4, 4, 3, 5, 5, 6};
    std::vector<std::vector<std::size_t>> const pigeon_ks {
        k_bl, k_ch, k_hp, k_ws};

    std::vector<std::size_t> const l_bl {11, 10, 8, 8, 10, 14, 11};
    std::vector<std::size_t> const l_ch {7, 11, 10, 9, 11, 9, 10, 10};
    std::vector<std::size_t> const l_hp {9, 12, 6, 12, 10, 11, 11, 9};
    std::vector<std::size_t> const l_ws {10, 9, 11, 13, 11, 11, 11, 12};
    std::vector<std::vector<std::size_t>> const pigeon_ls {
        k_bl, k_ch, k_hp, k_ws};

    for (unsigned site_id = 0U; site_id < site_paths.size(); ++site_id) {
        auto const& site = site_paths[site_id];
        auto const n_pigeons = sites[site_id].size();

        // read curves
        std::vector<Curves> pigeon_curves;
        pigeon_curves.reserve(n_pigeons);
        for (unsigned pid = 0U; pid < n_pigeons; ++pid) {
            auto const path = "data/Data_for_Mann_et_al_RSBL/" + site + "/"
                + sites[site_id][pid] + "/utm";
            Curves raw_curves = io::read_curves(path, 1);
            Curves curves;
            curves.reserve(raw_curves.size());
            for (auto const& curve: raw_curves)
                curves.emplace_back(curve.naive_l_simplification(200));
            pigeon_curves.push_back(std::move(curves));
        }

        // learning
        Curves centers;
        for (unsigned pid = 0U; pid < n_pigeons; ++pid) {
            auto const k = pigeon_ks[site_id][pid];
            auto const l = pigeon_ls[site_id][pid];
            auto clustering = clustering::computeCenterClustering(
                pigeon_curves[pid], k, l, true, true, cluster_alg, center_alg,
                CurveSimpMatrix({}, {}, df::frechet), init_dist, dist,
                df::frechet_lt, 1);
            for (auto& cluster: clustering)
                centers.push_back(std::move(cluster.center_curve));

            // for simplicity, fill up centers using the first such that there
            // are k in total.
            if (clustering.empty())
                ERROR("The clustering has to contain at least one center.");

            auto representative = centers.back();
            while (centers.size() % k != 0)
                centers.push_back(representative);
        }

        // testing on all trajectories
        std::size_t correct = 0;
        std::size_t wrong = 0;
        auto start = hrc::now();
        for (unsigned pid = 0U; pid < n_pigeons; ++pid) {
            for (auto const& curve: pigeon_curves[pid]) {
                // classify curve
                distance_t min_dist = std::numeric_limits<distance_t>::max();
                std::size_t min_id = 0;
                for (CurveID cid = 0; cid < centers.size(); ++cid) {
                    auto const& center = centers[cid];
                    distance_t d = dist(curve, center);
                    if (d < min_dist) {
                        min_dist = d;
                        min_id = cid;
                    }
                }

                // check if classification is correct
                auto const k = pigeon_ks[site_id][pid];
                if (pid == min_id / k)
                    ++correct;
                else {
                    ++wrong;
                    std::cout << "Pigeon " << pid << ":";
                    for (auto const& p: curve.get_points())
                        std::cout << " " << p.x << " " << p.y;
                    std::cout << "\nCenter:";
                    for (auto const& p: centers[min_id].get_points())
                        std::cout << " " << p.x << " " << p.y;
                    std::cout << "\n" << std::endl;
                }
            }
        }
        auto time = std::chrono::duration_cast<ms>(hrc::now() - start);

        std::cout << "Site: " << site_id << "\nTime: " << time.count()
                  << " ms\nCorrect classifications: " << correct
                  << "\nWrong classifications: " << wrong << "\nAccuracy: "
                  << static_cast<double>(correct) / (correct + wrong)
                  << std::endl;
    }
}
