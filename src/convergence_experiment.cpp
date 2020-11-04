#include "convergence_experiment.h"

#include <string>
#include "clustering/center_clustering_algs.h"
#include "distance_functions.h"
#include "utils/io.h"

namespace {
    void test_convergence(std::string const& base_path, std::size_t k,
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

void experiments::test_convergence() {
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
    for (std::size_t site_id = 0; site_id < site_directories.size(); ++site_id) {
        auto const& site_dir = site_directories[site_id];
        for (std::size_t pigeon_id = 0; pigeon_id < pigeon_directories[site_id].size(); ++pigeon_id) {
            auto const path = site_dir + pigeon_directories[site_id][pigeon_id];
            auto const k = pigeon_ks[site_id][pigeon_id];
            auto const l = pigeon_ls[site_id][pigeon_id];
            std::cout << pigeon_directories[site_id][pigeon_id] << ": ";
            test_convergence(path, k, l);
        }
    }
}
