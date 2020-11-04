#ifndef CENTER_CLUSTERING_ALGS_H
#define CENTER_CLUSTERING_ALGS_H

#include <cstddef>
#include <functional>
#include <string>
#include "clustering/center_algs.h"
#include "clustering/init_clustering_algs.h"
#include "Curve.h"

namespace detail {
    template<typename... Args>
    Clustering computeCenterClusteringRound(Curves const& curves,
            std::size_t k, std::size_t l, bool fix_start, bool fix_end,
            clustering::ClusterAlg cluster_alg, clustering::CenterAlg center_alg,
            DistanceMatrix<distance_t> const& dist_matrix,
            std::function<distance_t(Curve const&, Curve const&)> const&
                init_dist,
            std::function<distance_t(Curve const&, Curve const&)> const& dist,
            std::function<bool(Curve const&, Curve const&, distance_t)> const&
                lt_simp, Args&&... args) {
        auto clustering = clustering::computeClustering(curves, k, l,
            cluster_alg, dist_matrix, init_dist, lt_simp, false);

        // iterate as long as there are new centers
        unsigned count = 1;
        unsigned const max_count = 20;
        while (count <= max_count) {
            bool res = clustering::computeCenters(curves, clustering,
                center_alg, fix_start, fix_end, dist,
                std::forward<Args>(args)...);
            if (!res)
                break;
            clustering::updateClustering(clustering, curves, dist);
            ++count;
        }
        return clustering;
    }
}

namespace clustering {
    template<typename... Args>
    Clustering computeCenterClustering(Curves const& curves,
            std::size_t k, std::size_t l, bool fix_start, bool fix_end,
            ClusterAlg cluster_alg, CenterAlg center_alg,
            DistanceMatrix<distance_t> const& dist_matrix,
            std::function<distance_t(Curve const&, Curve const&)> const&
                init_dist,
            std::function<distance_t(Curve const&, Curve const&)> const& dist,
            std::function<bool(Curve const&, Curve const&, distance_t)> const&
                lt_simp,
            unsigned max_rounds, Args&&... args) {
        Clustering min_clustering;
        distance_t min_cost = std::numeric_limits<distance_t>::max();

        // To accommodate for randomness in initial clustering and updates, we
        // try the entire pipeline several times and pick the best result.
        for (unsigned round = 0; round < max_rounds; ++round) {
            auto clustering = detail::computeCenterClusteringRound(curves, k, l,
                fix_start, fix_end, cluster_alg, center_alg, dist_matrix,
                init_dist, dist, lt_simp, std::forward<Args>(args)...);
            distance_t cost_sum = 0.0;
            for (auto const& cluster: clustering)
                cost_sum += cluster.cost;

            if (cost_sum < min_cost) {
                min_clustering = std::move(clustering);
                min_cost = cost_sum;
            }
        }

        // Remove empty clusters.
        for (ClusterID cid = 0; cid < min_clustering.size(); ++cid) {
            if (min_clustering[cid].curve_ids.empty()) {
                std::swap(min_clustering[cid], min_clustering.back());
                min_clustering.pop_back();
            }
        }
        return min_clustering;
    }
}
#endif
