#pragma once

#include "Curve.h"
#include "distance_functions.h"
#include "clustering/clustering_algs.h"
#include "clustering/center_algs.h"
#include "clustering/center_clustering_algs.h"
#include "io.h"
#include "synthetic_curves.h"
#include "simplification/imaiiri.h"

#include <algorithm>
#include <chrono>
#include <string>
#include <unordered_map>
#include <vector>

// timing
using hrc = std::chrono::high_resolution_clock;
using time_point = hrc::time_point;
using ns = std::chrono::nanoseconds;



std::vector<Curve> sample(std::vector<Curve> curves, unsigned period) {

    std::vector<Curve> samples = std::vector<Curve>();

    for (std::size_t i = 0; i < curves.size(); ++i) {
        if (i % period == 0) {
            samples.push_back(curves[i]);
        }
    }

    return samples;
}



void characterClassification()
{
	int header_size = 0;
	ClusterAlg cluster_alg = ClusterAlg::Gonzalez;
	CenterAlg center_alg = CenterAlg::fCenter;
    
	int const k = 6;
	int const l = 12;
	int const k_xv = 5; // the k of the cross-validation

	auto curves = sample(io::read_curves("data/characters/data"), 15);
    std::cout << "loaded curves...\n";
	std::random_shuffle(curves.begin(), curves.end());
	FrechetLight frechet_light;

	// put curves into own vectors
	std::unordered_map<char, Curves> char_to_curves;
	for (auto const& curve: curves) {
		char_to_curves[curve.name().front()].push_back(curve);
	}
	std::vector<std::vector<Curves>> curves_by_char;
	for (auto const& pair: char_to_curves) {
		curves_by_char.emplace_back(k_xv);
		auto& curves = pair.second;
		for (CurveID i = 0; i < curves.size(); ++i) {
			curves_by_char.back()[i%k_xv].push_back(std::move(curves[i]));
		}
	}

	int correct = 0;
	int wrong = 0;
	std::size_t time = 0;
	for (int i = 0; i < k_xv; ++i) {
		Curves centers;
		for (auto const& curves_slots: curves_by_char) {
			Curves training_curves;
			for (int j = 0; j < k_xv; ++j) {
				if (i == j) { continue; }
				training_curves.insert(
					training_curves.end(), curves_slots[j].begin(), curves_slots[j].end());
			}
			auto clustering = computeCenterClustering(
				training_curves, k, l, cluster_alg, center_alg, frechet, "", 1);
			for (auto const& cluster: clustering) {
				centers.push_back(std::move(cluster.center_curve));
			}

			// for simplicity, fill up centers using the first such that there are k in total.
			if (clustering.empty()) { ERROR("The clustering has to contain at least one center."); }
			auto representative = centers.back();
			while (centers.size()%k != 0) { centers.push_back(representative); }
		}

		// cross-validate on remaining data
		auto start = hrc::now();
		for (std::size_t char_id = 0; char_id < curves_by_char.size(); ++char_id) {
			auto const& curves = curves_by_char[char_id][i];
			for (auto const& curve: curves) {
				// classify curve
				distance_t min_dist = std::numeric_limits<distance_t>::max();
				std::size_t min_id = 0;
				for (CurveID center_id = 0; center_id < centers.size(); ++center_id) {
					auto const& center = centers[center_id];
					// if (frechet_light.lessThanWithFilters(min_dist, curve, center)) {
                    distance_t dist = frechet(curve, center);
					if (dist < min_dist) {
                    	// auto dist = frechet_light.calcDistance(curve, center);
						min_dist = dist;
						min_id = center_id;
					}
				}

				// check if classification is correct
				if (char_id == min_id/k) ++correct; else ++wrong;
			}
		}
		time += std::chrono::duration_cast<ns>(hrc::now()-start).count();
	}

	std::cout << "Time: " << time/1000000. << " ms\n";
	std::cout << "Correct classifications: " << correct << "\n";
	std::cout << "Wrong classifications: " << wrong << "\n";
	std::cout << "Success rate: " << float(correct)/(correct+wrong) << "\n";
}

