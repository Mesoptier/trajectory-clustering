#include "convergence_experiment.h"

#include <string>
#include "clustering/center_clustering_algs.h"
#include "distance_functions.h"
#include "utils/io.h"

namespace {
    void testConvergence(std::string const& base_path, std::size_t k,
            std::size_t l) {
        auto curves = io::read_curves(base_path, 1);

        auto cluster_alg = clustering::ClusterAlg::PAM;
        auto center_alg = clustering::CenterAlg::cdba;

        auto clustering = clustering::computeCenterClustering(curves, k, l,
            true, true, cluster_alg, center_alg,
            CurveSimpMatrix({}, {}, df::frechet), df::integral_frechet,
            df::integral_frechet, df::frechet_lt, 1);

        // calculate max_filename_id
        std::vector<std::pair<std::string, ClusterID>> filename_id_pairs;
        for (std::size_t cluster_id = 0; cluster_id < clustering.size(); ++cluster_id) {
            auto const& cluster = clustering[cluster_id];
            for (auto curve_id: cluster.curve_ids) {
                filename_id_pairs.emplace_back(curves[curve_id].filename, cluster_id);
            }
        }
        std::sort(filename_id_pairs.begin(), filename_id_pairs.end());
        std::string max_filename = filename_id_pairs.back().first;
        ClusterID max_filename_id = filename_id_pairs.back().second;

        // calculate max_size_id
        auto size_comp = [&](Cluster const& c1, Cluster const& c2) {
            return c1.curve_ids.size() < c2.curve_ids.size();
        };
        auto it = std::max_element(clustering.begin(), clustering.end(), size_comp);
        ClusterID max_size_id =
            static_cast<std::size_t>(std::distance(clustering.begin(), it));

        if (max_filename_id == max_size_id) {
            // Additionally calculate suffix length
            long i = static_cast<long>(filename_id_pairs.size()) - 1;
            while (i >= 0) {
                if (filename_id_pairs[static_cast<std::size_t>(i)].second != max_size_id) { break; }
                --i;
            }
            auto suffix_length = static_cast<long>(filename_id_pairs.size()) - 1 - i;
            std::cout << "Yes! Suffix length: " << suffix_length << "/" << curves.size() << "\n";
        }
        else {
            std::cout << "No!" << "\n";
        }
    }
}

void experiments::testConvergence() {
    for (std::size_t site_id = 0; site_id < site_directories.size(); ++site_id) {
        auto const& site_dir = site_directories[site_id];
        for (std::size_t pigeon_id = 0; pigeon_id < pigeon_directories[site_id].size(); ++pigeon_id) {
            auto const path = site_dir + pigeon_directories[site_id][pigeon_id];
            auto const k = pigeon_ks[site_id][pigeon_id];
            auto const l = pigeon_ls[site_id][pigeon_id];
            std::cout << pigeon_directories[site_id][pigeon_id] << ": ";
            testConvergence(path, k, l);
        }
    }
}
