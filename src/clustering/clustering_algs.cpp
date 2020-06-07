#include "clustering/clustering_algs.h"

#include <algorithm>
#include <limits>
#include <numeric>
#include <unordered_map>

#include "simplification/agarwal.h"
#include "utils/defs.h"
#include "utils/matrix.h"
#include "utils/random.h"
#include "utils/union_find.h"

namespace {
    // TODO: Computes all distances, not only one per pair.
    template <typename Comp>
    Clustering linkage(Curves const& curves,
            std::size_t k, std::size_t l, Comp comp,
            std::function<distance_t(Curve const&, Curve const&)> const& dist,
            bool naive_simplification) {
        // compute all pairwise distances
        Matrix<distance_t> dist_matrix(curves.size(), curves.size());
        for (CurveID cid1 = 0; cid1 < curves.size(); ++cid1) {
            for (CurveID cid2(cid1); cid2 < curves.size(); ++cid2) {
                if (cid1 == cid2)
                    dist_matrix(cid1, cid2) = 0;
                else {
                    auto d = dist(curves[cid1], curves[cid2]);
                    dist_matrix(cid1, cid2) = d;
                    dist_matrix(cid2, cid1) = d;
                }
            }
        }

        // create initial clusters
        CurveIDs base_set(curves.size());
        std::iota(base_set.begin(), base_set.end(), 0);
        UnionFind<CurveID> union_find(base_set);

        // merge clusters until there are exactly k
        while (union_find.getRoots().size() > k) {
            // find two clusters to merge
            distance_t min_dist = std::numeric_limits<distance_t>::max();
            CurveID min_id1, min_id2;
            for (auto const& curve_id1: union_find.getRoots()) {
                for (auto const& curve_id2: union_find.getRoots()) {
                    if (curve_id1 == curve_id2)
                        continue;

                    auto new_dist = dist_matrix(curve_id1, curve_id2);
                    if (new_dist < min_dist) {
                        min_dist = new_dist;
                        min_id1 = curve_id1;
                        min_id2 = curve_id2;
                    }
                }
            }

            // merge clusters and adapt distances
            min_id1 = union_find.findRoot(min_id1);
            min_id2 = union_find.findRoot(min_id2);
            auto new_root_id = union_find.uniteSets(min_id1, min_id2);
            auto new_child_id = (new_root_id == min_id1 ? min_id2 : min_id1);

            for (auto const& curve_id: union_find.getRoots()) {
                auto root_dist = dist_matrix(new_root_id, curve_id);
                auto child_dist = dist_matrix(new_child_id, curve_id);
                dist_matrix(new_root_id, curve_id) = comp(root_dist, child_dist);
            }
        }

        // construct the result
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
                simplification::greedy::simplify(curves[curve_id], l,
                    [dist](const Curve& a, const Curve& b, distance_t t) {
                        return dist(a, b) <= t;
                });
        }

        return result;
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
        std::function<distance_t(Curve const&, Curve const&)> const& dist,
        std::string const& dist_matrix, bool naive_simplification) {
    switch (cluster_alg) {
    case ClusterAlg::SingleLinkage:
        return singleLinkage(curves, k, l, dist, naive_simplification);
    case ClusterAlg::CompleteLinkage:
        return completeLinkage(curves, k, l, dist, naive_simplification);
    case ClusterAlg::Gonzalez:
        return gonzalez(curves, k, l, dist, naive_simplification);
    case ClusterAlg::PAM:
        return pam_with_simplifications(curves, k, l, dist, dist_matrix);
    }
    ERROR("No matching cluster_alg enum passed.");
}

Clustering clustering::singleLinkage(Curves const& curves,
        std::size_t k, std::size_t l,
        std::function<distance_t(Curve const&, Curve const&)> const& dist,
        bool naive_simplification) {
    auto min = [](distance_t a, distance_t b) {
        return std::min(a,b);
    };
    return linkage(curves, k, l, min, dist, naive_simplification);
}

Clustering clustering::completeLinkage(Curves const& curves,
        std::size_t k, std::size_t l,
        std::function<distance_t(Curve const&, Curve const&)> const& dist,
        bool naive_simplification) {
    auto max = [](distance_t a, distance_t b) {
        return std::max(a, b);
    };
    return linkage(curves, k, l, max, dist, naive_simplification);
}

Clustering clustering::gonzalez(Curves const& curves,
        std::size_t k, std::size_t l,
        std::function<distance_t(Curve const&, Curve const&)> const& dist,
        bool naive_simplification) {
    Clustering result;

    auto max_dist = std::numeric_limits<distance_t>::max();
    std::vector<distance_t> distances_to_center(curves.size(), max_dist);
    ClusterIDs closest_center(curves.size());

    Random random;
    CurveID center_id =
        random.getUniformInt<Curves::size_type>(0, curves.size() - 1);

    // add as center and update closest distances to center
    auto center_curve = naive_simplification ?
        curves[center_id].naive_l_simplification(l) :
        simplification::greedy::simplify(curves[center_id], l,
            [dist](const Curve& a, const Curve& b, distance_t t) {
                return dist(a, b) <= t;
        });
    
    result.push_back({{}, center_curve});
    for (CurveID curve_id = 0; curve_id < curves.size(); ++curve_id) {
        auto& current_dist = distances_to_center[curve_id];
        auto new_dist = dist(center_curve, curves[curve_id]);
        if (new_dist < current_dist) {
            current_dist = new_dist;
            closest_center[curve_id] = result.size() - 1;
        }
    }

    while (result.size() < k) {
        auto center_it = std::max_element(distances_to_center.begin(),
                                          distances_to_center.end());
        auto cid = static_cast<std::size_t>(
            std::distance(distances_to_center.begin(), center_it));
        Curve cent_curve = naive_simplification ? 
            curves[cid].naive_l_simplification(l) : 
            simplification::greedy::simplify(curves[cid], l,
                [dist](const Curve& a, const Curve& b, distance_t t) {
                    return dist(a, b) <= t;
            });

        result.push_back({{}, cent_curve});
        for (CurveID curve_id = 0; curve_id < curves.size(); ++curve_id) {
            auto& current_dist = distances_to_center[curve_id];

            auto new_dist = dist(cent_curve, curves[curve_id]);
            if (new_dist < current_dist) {
                current_dist = new_dist;
                closest_center[curve_id] = result.size() - 1;
            }
        }
    }

    for (CurveID curve_id = 0; curve_id < curves.size(); ++curve_id) {
        auto cluster_id = closest_center[curve_id];
        result[cluster_id].curve_ids.push_back(curve_id);
    }

    return result;
}

void clustering::updateClustering(Curves const& curves, Clustering& clustering,
        std::function<distance_t(Curve const&, Curve const&)> const& dist) {
    // clear clusters
    for (auto& cluster: clustering)
        cluster.curve_ids.clear();

    // compute new clusters
    for (CurveID curve_id = 0; curve_id < curves.size(); ++curve_id) {
        distance_t min_dist = std::numeric_limits<distance_t>::max();
        ClusterID min_cluster_id;
        for (ClusterID cid = 0; cid < clustering.size(); ++cid) {
            auto const& center_curve = clustering[cid].center_curve;

            auto new_dist = dist(curves[cid], center_curve);
            if (new_dist < min_dist) {
                min_dist = new_dist;
                min_cluster_id = cid;
            }
        }
        clustering[min_cluster_id].curve_ids.push_back(curve_id);
    }
}

distance_t clustering::calcDiameter(Curves const& curves,
        std::function<distance_t(Curve const&, Curve const&)> const& dist) {
    distance_t max_distance = 0.0;
    for (CurveID cid1 = 0; cid1 < curves.size(); ++cid1)
        for (CurveID cid2(cid1 + 1); cid2 < curves.size(); ++cid2)
            max_distance = std::max(max_distance,
                                    dist(curves[cid1], curves[cid2]));

    return max_distance;
}

Clustering clustering::pam_with_simplifications(Curves const& curves,
        std::size_t k, std::size_t l,
        std::function<distance_t(Curve const&, Curve const&)> const& dist,
        std::string const& matrix_file_name, bool naive_simplification) {
    Curves simplifications;
    auto lt = [dist](Curve const& a, Curve const& b, distance_t t) {
        return dist(a, b) <= t;
    };

    for (auto const& curve: curves)
        simplifications.emplace_back(naive_simplification ?
            curve.naive_l_simplification(l) :
            simplification::greedy::simplify(curve, l, lt));

    CurveSimpMatrix distance_matrix = matrix_file_name == "" ? 
        CurveSimpMatrix(curves, simplifications, dist) :
        CurveSimpMatrix(matrix_file_name);

    std::vector<std::size_t> initial_centers =
        clustering::pam::compute(curves.size(), k, distance_matrix);

    Clustering clustering;
    for (auto const& center_id: initial_centers) {
        DEBUG(curves[center_id].name());
        clustering.push_back({{}, simplifications[center_id], 0});
    }
    updateClustering(curves, clustering, dist);
    return clustering;
}

Clustering clustering::pam_with_centering(Curves const& curves,
        std::size_t k, std::size_t l,
        std::function<distance_t(Curve const&, Curve const&)> const& dist,
        std::string const& matrix_file_name, bool naive_simplification) {
    Clustering clustering = pam_with_simplifications(curves, k, l, dist,
        matrix_file_name, naive_simplification);

    for (auto& cluster: clustering) {
        cluster.cost = calcC2CDist(curves, cluster.center_curve,
            cluster.curve_ids, C2CDist::Median, dist);
    }

    calcFSACenters(curves, clustering, l, dist, C2CDist::Median,
        CenterCurveUpdateMethod::frechetMean);
    updateClustering(curves, clustering, dist);
    return clustering;
}

distance_t clustering::kMedianCost(Curves const& curves,
        Clustering const& clustering,
        std::function<distance_t(Curve const&, Curve const&)> const& dist) {
    distance_t cost = 0.0;
    for (auto const& cluster: clustering)
        cost += calcC2CDist(curves, cluster.center_curve, cluster.curve_ids,
            C2CDist::Median, dist);
    return cost;
}
