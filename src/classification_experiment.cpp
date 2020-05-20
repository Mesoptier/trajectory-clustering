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

Curves classification::sample(const Curves& curves, unsigned period) {
    Curves samples;
    for (std::size_t i = 0; i < curves.size(); i += period)
        samples.push_back(curves[i]);
    return samples;
}

void classification::characterClassification(
        clustering::ClusterAlg cluster_alg, clustering::CenterAlg center_alg) {
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
    for (auto const& curve: curves) {
        char_to_curves[curve.name().front()].push_back(curve);
    }

    // Upper level: characters, second level: CV folds
    std::vector<std::vector<Curves>> curves_by_char;
    for (auto const& pair: char_to_curves) {
        curves_by_char.emplace_back(k_xv);
        auto& cvs = pair.second;
        for (CurveID i = 0; i < cvs.size(); ++i) {
            curves_by_char.back()[i%k_xv].push_back(std::move(cvs[i]));
        }
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
            auto clustering = computeCenterClustering(training_curves, k, l,
                cluster_alg, center_alg, df::frechet, df::frechet, "", 1);
            for (auto const& cluster: clustering)
                centers.push_back(std::move(cluster.center_curve));

            // for simplicity, fill up centers using the first such that there
            // are k in total.
            if (clustering.empty()) {
                ERROR("The clustering has to contain at least one center.");
            }
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
                    distance_t dist = df::frechet(curve, center);
                    if (dist < min_dist) {
                        min_dist = dist;
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

    std::cout << "Time: " << std::chrono::duration_cast<ms>(time).count() << " ms\n";
    std::cout << "Correct classifications: " << correct << "\n";
    std::cout << "Wrong classifications: " << wrong << "\n";
    std::cout << "Accuracy: " << static_cast<double>(correct) /
        (correct + wrong) << "\n";
}
