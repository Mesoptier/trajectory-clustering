#include "clustering/center_clustering_algs.h"

#include <utility>

namespace {
    Clustering computeCenterClusteringRound(Curves const& curves,
            std::size_t k, std::size_t l,
            clustering::ClusterAlg cluster_alg, clustering::CenterAlg center_alg,
            std::function<distance_t(Curve const&, Curve const&)> const& dist,
            std::string const& dist_matrix) {
        auto clustering = computeClustering(curves, k, l, cluster_alg, dist,
            dist_matrix, true);
        clustering::updateClustering(curves, clustering, dist);

        // iterate as long as there are new centers
        unsigned count = 1;
        unsigned const max_count = 10;
        while (count <= max_count && clustering::computerCenters(
                curves, clustering, l, center_alg, dist)) {
            clustering::updateClustering(curves, clustering, dist);
            ++count;
        }
        return clustering;
    }
}

Clustering clustering::computeCenterClustering(Curves const& curves,
        std::size_t k, std::size_t l,
        ClusterAlg cluster_alg, CenterAlg center_alg,
        std::function<distance_t(Curve const&, Curve const&)> const& dist,
        std::string const& dist_matrix, unsigned max_rounds) {
    Clustering min_clustering;
    distance_t min_cost = std::numeric_limits<distance_t>::max();

    for (unsigned round = 0; round < max_rounds; ++round) {
        auto clustering = computeCenterClusteringRound(curves, k, l,
            cluster_alg, center_alg, dist, dist_matrix);
        distance_t cost_sum = 0.0;
        for (auto const& cluster: clustering)
            cost_sum += cluster.cost;

        if (cost_sum < min_cost) {
            min_clustering = std::move(clustering);
            min_cost = cost_sum;
        }
    }

    // remove empty clusters
    for (ClusterID cid = 0; cid < min_clustering.size(); ++cid) {
        if (min_clustering[cid].curve_ids.empty()) {
            std::swap(min_clustering[cid], min_clustering.back());
            min_clustering.pop_back();
        }
    }
    return min_clustering;
}
