#include "clustering/init_clustering_algs.h"

#include <algorithm>
#include <limits>
#include <numeric>
#include <unordered_map>

#include "geom.h"
#include "clustering/pam.h"
#include "simplification/agarwal.h"
#include "utils/defs.h"
#include "utils/random.h"
#include "utils/union_find.h"

namespace {
    /**
     * \brief Linkage clustering framework.
     *
     * Starts with one curve per cluster and merges them until we get k
     * clusters. Single linkage vs complete linkage: different function for
     * comp, so different distances between clusters after merging.
     * \param curves The curves to cluster.
     * \param k The number of clusters.
     * \param l The complexity of centers.
     * \param comp The linkage function.
     * \param dist_matrix The distance matrix between all non-simplified curves.
     * \param lt_simp The distance function to use for simplifying centers.
     * \param naive_simplification Whether to use naive simplification or the
     * greedy approach with lt_simp.
     * \return The clustering.
     */
    template <typename Comp>
    Clustering linkage(Curves const& curves, std::size_t k, std::size_t l,
            Comp comp, SymmetricMatrix dist_matrix,
            std::function<bool(Curve const&, Curve const&, distance_t)> const&
            lt_simp, bool naive_simplification) {
        // Create initial clusters.
        CurveIDs base_set(curves.size());
        std::iota(base_set.begin(), base_set.end(), 0);
        UnionFind<CurveID> union_find(base_set);

        // Merge clusters until there are exactly k.
        while (union_find.getRoots().size() > k) {
            // Find two clusters to merge.
            distance_t min_dist = std::numeric_limits<distance_t>::max();
            CurveID min_id1, min_id2;
            for (auto const& curve_id1: union_find.getRoots()) {
                for (auto const& curve_id2: union_find.getRoots()) {
                    if (curve_id1 == curve_id2)
                        continue;

                    auto new_dist = dist_matrix.at(curve_id1, curve_id2);
                    if (new_dist < min_dist) {
                        min_dist = new_dist;
                        min_id1 = curve_id1;
                        min_id2 = curve_id2;
                    }
                }
            }

            // Merge clusters and adapt distances.
            min_id1 = union_find.findRoot(min_id1);
            min_id2 = union_find.findRoot(min_id2);
            auto new_root_id = union_find.uniteSets(min_id1, min_id2);
            auto new_child_id = (new_root_id == min_id1 ? min_id2 : min_id1);

            for (auto const& curve_id: union_find.getRoots()) {
                auto root_dist = dist_matrix.at(new_root_id, curve_id);
                auto child_dist = dist_matrix.at(new_child_id, curve_id);
                dist_matrix.at(new_root_id, curve_id) = 
                    comp(root_dist, child_dist);
            }
        }

        // Construct the result.
        Clustering result(union_find.getRoots().size());

        std::unordered_map<CurveID, ClusterID> to_cluster_id;
        ClusterID cluster_id = 0;
        for (auto const& curve_id: union_find.getRoots()) {
            to_cluster_id[curve_id] = cluster_id;
            ++cluster_id;
        }
        for (CurveID curve_id = 0; curve_id < curves.size(); ++curve_id) {
            auto cl_id = to_cluster_id[union_find.findRoot(curve_id)];
            result[cl_id].curve_ids.push_back(curve_id);
        }

        // We just take the root curves as centers. They don't have any special
        // meaning, but at least we supply some centers.
        for (auto const& curve_id: union_find.getRoots()) {
            auto cl_id = to_cluster_id[curve_id];
            result[cl_id].center_curve = naive_simplification ?
                curves[curve_id].naive_l_simplification(l) :
                simplification::greedy::simplify(curves[curve_id], l, lt_simp);
            result[cl_id].center_id = curve_id;
        }
        return result;
    }

    /**
     * \brief Update the distances to the closest center after adding a new
     * center, as used in Gonzales' algorithm.
     *
     * If using the distance matrix, make sure the simplifications are computed
     * in the same way as implied by naive_simplification and lt_simp.
     * \param curves The curves being clustered.
     * \param l The complexity of centers.
     * \param dist_matrix The distance matrix between curves and their
     * simplifications. May be empty.
     * \param use_mtrx Whether to use the dist_matrix or dist.
     * \param dist The distance function. Not used if use_mtrx is true.
     * \param lt_simp The distance function for the greedy simplification.
     * \param naive_simplification Whether to use naive or greedy simplification.
     * \param cid The new center id.
     * \param closest_center The vector indicating the closest center id for
     * each curve. Updated when adding a new center.
     * \param distances_to_center The vector indicating the distance to the
     * closest center for each curve. Updated when adding a new center.
     * \param result The resulting clustering.
     */
    void add_center_curve_update_distances(Curves const& curves,
            std::size_t l, CurveSimpMatrix const& dist_matrix, bool use_mtrx,
            std::function<distance_t(Curve const&, Curve const&)> const& dist,
            std::function<bool(Curve const&, Curve const&, distance_t)> const&
            lt_simp, bool naive_simplification,
            CurveID cid, ClusterIDs& closest_center,
            std::vector<distance_t>& distances_to_center, Clustering& result) {
        auto center_curve = naive_simplification ?
            curves[cid].naive_l_simplification(l) : 
            simplification::greedy::simplify(curves[cid], l, lt_simp);

        result.push_back({{}, center_curve, cid});
        distances_to_center[cid] = 0.0;
        closest_center[cid] = result.size() - 1;

        for (CurveID curve_id = 0; curve_id < curves.size(); ++curve_id) {
            auto& current_dist = distances_to_center[curve_id];
            auto new_dist = use_mtrx ? dist_matrix.at(curve_id, cid)
                : dist(center_curve, curves[curve_id]);
            if (new_dist < current_dist) {
                current_dist = new_dist;
                closest_center[curve_id] = result.size() - 1;
            }
        }
    }

    /**
     * \brief Update the lists of curve ids in each cluster.
     * \param clustering The clustering with correct centers.
     * \param dist_matrix The distance matrix for curves and simplifications.
     */
    void updateClustering(Clustering& clustering,
            CurveSimpMatrix const& dist_matrix) {
        // Clear clusters.
        for (auto& cluster: clustering)
            cluster.curve_ids.clear();

        // Compute new clusters.
        for (CurveID curve_id = 0; curve_id < dist_matrix.size(); ++curve_id) {
            auto min_dist = std::numeric_limits<distance_t>::max();
            ClusterID min_cluster_id;
            for (ClusterID cid = 0; cid < clustering.size(); ++cid) {
                auto new_dist = dist_matrix.at(curve_id,
                    clustering[cid].center_id);
                if (new_dist < min_dist) {
                    min_dist = new_dist;
                    min_cluster_id = cid;
                }
            }
            clustering[min_cluster_id].curve_ids.push_back(curve_id);
        }
    }
} // end anonymous namespace

std::string clustering::toString(ClusterAlg cluster_alg) {
    switch(cluster_alg) {
    case ClusterAlg::SingleLinkage: return "SingleLinkage";
    case ClusterAlg::CompleteLinkage: return "CompleteLinkage";
    case ClusterAlg::Gonzalez: return "Gonzalez";
    case ClusterAlg::PAM: return "PAM";
    }
    ERROR("Unknown cluster_alg.");
}

Clustering clustering::computeClustering(Curves const& curves,
        std::size_t k, std::size_t l, ClusterAlg cluster_alg,
        DistanceMatrix<distance_t> const& dist_matrix,
        std::function<distance_t(Curve const&, Curve const&)> const& dist,
        std::function<bool(Curve const&, Curve const&, distance_t)> const& lt_simp,
        bool naive_simplification) {
    switch (cluster_alg) {
    case ClusterAlg::SingleLinkage:
        return singleLinkage(curves, k, l,
            dynamic_cast<SymmetricMatrix const&>(dist_matrix),
            lt_simp, naive_simplification);
    case ClusterAlg::CompleteLinkage:
        return completeLinkage(curves, k, l,
            dynamic_cast<SymmetricMatrix const&>(dist_matrix),
            lt_simp, naive_simplification);
    case ClusterAlg::Gonzalez:
        return gonzalez(curves, k, l,
            dynamic_cast<CurveSimpMatrix const&>(dist_matrix), dist,
            lt_simp, naive_simplification);
    case ClusterAlg::PAM:
        return pam_with_simplifications(curves, k, l,
            dynamic_cast<CurveSimpMatrix const&>(dist_matrix),
            lt_simp, naive_simplification);
    }
    ERROR("No matching cluster_alg enum passed.");
}

Clustering clustering::singleLinkage(Curves const& curves,
        std::size_t k, std::size_t l, SymmetricMatrix const& dist_matrix,
        std::function<bool(Curve const&, Curve const&, distance_t)> const& lt_simp,
        bool naive_simplification) {
    auto min = [](distance_t a, distance_t b) {
        return std::min(a, b);
    };
    return linkage(curves, k, l, min, dist_matrix, lt_simp, naive_simplification);
}

Clustering clustering::completeLinkage(Curves const& curves,
        std::size_t k, std::size_t l, SymmetricMatrix const& dist_matrix,
        std::function<bool(Curve const&, Curve const&, distance_t)> const& lt_simp,
        bool naive_simplification) {
    auto max = [](distance_t a, distance_t b) {
        return std::max(a, b);
    };
    return linkage(curves, k, l, max, dist_matrix, lt_simp, naive_simplification);
}

Clustering clustering::gonzalez(Curves const& curves,
        std::size_t k, std::size_t l, CurveSimpMatrix const& dist_matrix,
        std::function<distance_t(Curve const&, Curve const&)> const& dist,
        std::function<bool(Curve const&, Curve const&, distance_t)> const& lt_simp,
        bool naive_simplification) {
    Clustering result;
    std::vector<distance_t> distances_to_center(curves.size(),
        std::numeric_limits<distance_t>::max());
    ClusterIDs closest_center(curves.size());
    bool use_mtrx = dist_matrix.size() == curves.size();

    // 1. Pick a random center.
    Random random;
    CurveID cid = random.getUniformInt<Curves::size_type>(0, curves.size() - 1);

    // 2. Compute distances from all curves to the cluster center.
    add_center_curve_update_distances(curves, l, dist_matrix, use_mtrx, dist,
        lt_simp, naive_simplification, cid,
        closest_center, distances_to_center, result);

    while (result.size() < k) {
        // 3. Pick a curve that is furthest away from any current center as the
        // new center and update the distances and clusters.
        auto center_it = std::max_element(distances_to_center.begin(),
                                          distances_to_center.end());
        cid = static_cast<std::size_t>(
            std::distance(distances_to_center.begin(), center_it));
        add_center_curve_update_distances(curves, l, dist_matrix, use_mtrx,
            dist, lt_simp, naive_simplification, cid,
            closest_center, distances_to_center, result);
    }

    // 4. Compute the clusters associated with the final centers.
    for (CurveID curve_id = 0; curve_id < curves.size(); ++curve_id) {
        auto cluster_id = closest_center[curve_id];
        result[cluster_id].curve_ids.push_back(curve_id);
    }
    return result;
}

Clustering clustering::pam_with_simplifications(Curves const& curves,
        std::size_t k, std::size_t l, CurveSimpMatrix const& dist_matrix,
        std::function<bool(Curve const&, Curve const&, distance_t)> const& lt_simp,
        bool naive_simplification) {
    std::vector<std::size_t> initial_centers =
        clustering::pam::compute(curves.size(), k, dist_matrix);

    // Convert the result into the required format and compute simplifications
    // of the center curves.
    Clustering clustering;
    for (auto const& center_id: initial_centers) {
        DEBUG(curves[center_id].name());
        auto simplified_center = naive_simplification ?
            curves[center_id].naive_l_simplification(l) :
            simplification::greedy::simplify(curves[center_id], l, lt_simp);
        clustering.push_back({{}, simplified_center, center_id, 0});
    }
    ::updateClustering(clustering, dist_matrix);
    return clustering;
}

void clustering::updateClustering(Clustering& clustering, Curves const& curves,
        std::function<distance_t(Curve const&, Curve const&)> const& dist) {
    // Clear clusters.
    for (auto& cluster: clustering)
        cluster.curve_ids.clear();

    // Compute new clusters.
    for (CurveID curve_id = 0; curve_id < curves.size(); ++curve_id) {
        auto min_dist = std::numeric_limits<distance_t>::max();
        ClusterID min_cluster_id;
        for (ClusterID cid = 0; cid < clustering.size(); ++cid) {
            auto new_dist = dist(curves[curve_id], clustering[cid].center_curve);
            if (new_dist < min_dist) {
                min_dist = new_dist;
                min_cluster_id = cid;
            }
        }
        clustering[min_cluster_id].curve_ids.push_back(curve_id);
    }
}

// distance_t clustering::kMedianCost(Curves const& curves,
//         Clustering const& clustering,
//         std::function<distance_t(Curve const&, Curve const&)> const& dist) {
//     distance_t cost = 0.0;
//     for (auto const& cluster: clustering)
//         cost += calcC2CDist(curves, cluster.center_curve, cluster.curve_ids,
//             C2CDist::Median, dist);
//     return cost;
// }
