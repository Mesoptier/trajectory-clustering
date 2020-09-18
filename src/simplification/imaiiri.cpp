#include "simplification/imaiiri.h"

#include <algorithm>
#include <limits>
#include <queue>
#include <unordered_set>
#include <vector>
#include "utils/SymmetricMatrix.h"

namespace {
    /**
     * \brief Dynamic program cell.
     */
    struct entry {
        distance_t cost;
        PointID prev;
        entry(): cost(std::numeric_limits<distance_t>::infinity()), prev(0) {}
    };

    /**
     * \brief Aggregate two values: either add up or return the maximum.
     * \param a The first value.
     * \param b The second value.
     * \param max Whether to use max or sum.
     * \return a + b or max(a, b).
     */
    inline distance_t aggr(distance_t a, distance_t b, bool max) {
        if (max)
            return a >= b ? a : b;
        return a + b;
    }

    /**
     * \brief Get the area of the bounding box of a curve to use as upper bound
     * for the threshold.
     * \param c The curve.
     * \return The area of the axis-aligned bounding box fitted to c.
     */
    distance_t bb_area(Curve const& c) {
        auto const& extreme = c.getExtremePoints();
        auto a = extreme.max_x - extreme.min_x;
        auto b = extreme.max_y - extreme.min_y;
        return a * b;
    }
}

std::pair<distance_t, Curve> simplification::imai_iri::simplify(Curve const& in,
        PointID const& ell,
        std::function<distance_t(Curve const&, Curve const&)> const& dist,
        bool max) {
    assert(in.size() > 1);
    assert(ell > 0);
    SymmetricMatrix distances(in.size());
    for (PointID i(0); i < in.size(); ++i) {
        for (PointID j(i); j < in.size(); ++j) {
            if (i == j)
                distances.at(i, j) = std::numeric_limits<distance_t>::infinity();
            else if (j == i + 1)
                distances.at(i, j) = 0.0;
            else
                distances.at(i, j) = dist(in.slice(i, j),
                                          Curve("", {in[i], in[j]}));
        }
    }

    // Q[k][i] = min cost to simplify in[0]...in[i] with k + 1 segments
    std::vector<std::vector<entry>> Q(ell,
        std::vector<entry>(in.size(), entry()));

    // Base case: first segment from p_0 to p_i, for all 1 <= i < n
    for (PointID i = 1; i < in.size(); ++i) {
        Q[0][i].cost = distances.at(0, i);
    }

    // Fill matrix row by row (increasing # segments)
    for (PointID k = 1; k < ell; ++k) {
        for (PointID i = k + 1; i < in.size(); ++i) {
            // Q[k][i] = min_{k <= j < i} (Q[k - 1][j] + d(segm, subcurve))
            // or: max(Q[k - 1][j], d(segm, subcurve))
            for (PointID j = k; j < i; ++j) {
                distance_t opt = aggr(Q[k - 1][j].cost, distances.at(j, i), max);
                if (opt < Q[k][i].cost) {
                    Q[k][i].cost = opt;
                    Q[k][i].prev = j;
                }
            }
        }
    }

    distance_t cost = std::numeric_limits<distance_t>::infinity();
    PointID best = ell;
    for (PointID k = 0; k < ell; ++k) {
        if (Q[k][in.size() - 1].cost < cost) {
            cost = Q[k][in.size() - 1].cost;
            best = k;
        }
    }
    assert(best < ell);
    Points simpl_points(best + 2);
    PointID i = in.size() - 1;
    simpl_points[best + 1] = in[i];
    while (best > 0) {
        auto tmp = in[Q[best][i].prev];
        simpl_points[best] = tmp;
        i = Q[best][i].prev;
        --best;
    }
    simpl_points[0] = in[0];
    return {cost, Curve("", simpl_points)};
}

Curve simplification::imai_iri::simplify(Curve const& in, PointID const& ell,
        std::function<bool(Curve const&, Curve const&, distance_t)> const&
        less_than) {
    // Note: if you know that the value of distance is the distance between
    // some two points / point and a segment, you can get the exact set of
    // possible values for t. As we do not make such assumption, we get the
    // general lower and upper bounds for t and do a numeric (binary) search.

    static constexpr distance_t epsilon = 1e-8;
    assert(in.size() > 1);
    assert(ell > 0);
    if (in.size() - 1 <= ell)
        return in;

    distance_t min = 0.0;
    distance_t max = bb_area(in);

    Curve simplified;
    while (max - min > epsilon) {
        auto split = (max + min) / 2.0;
        simplified = simplify(in, split, less_than);
        if (simplified.size() <= ell)
            max = split;
        else
            min = split;
    }

    return simplify(in, max, less_than);
}

Curve simplification::imai_iri::simplify(Curve const& in, distance_t delta,
        std::function<bool(Curve const&, Curve const&, distance_t)> const&
        less_than) {
    // Note: for Hausdorff or Fréchet distance we could cleverly check whether
    // disks are stabbed [Chin, Chan 1992] (in the right order for Fréchet
    // [Guibas et al. 1993]) and not call less_than, but we do not do that in
    // the interest of a more general solution.

    std::vector<std::unordered_set<PointID>> adj(in.size() - 1);

    // Shortcut graph
    for (PointID i = 0; i < in.size() - 1; ++i) {
        adj[i].emplace(i + 1);
        for (PointID j = i + 2; j < in.size(); ++j) {
            if (less_than(in.slice(i, j), Curve("", {in[i], in[j]}), delta))
                adj[i].emplace(j);
        }
    }

    // BFS
    std::vector<bool> visited(in.size(), false);
    std::vector<PointID> parent(in.size(), PointID::invalid_value);
    std::queue<PointID> bfs_q;
    bfs_q.emplace(0);
    visited[0] = true;
    while (!bfs_q.empty()) {
        auto v = bfs_q.front();
        bfs_q.pop();
        if (v == in.size() - 1)
            break;
        for (auto const& w: adj[v]) {
            if (!visited[w]) {
                visited[w] = true;
                parent[w] = v;
                bfs_q.emplace(w);
            }
        }
    }

    // Reconstruct the curve
    Points simpl;
    PointID i = in.size() - 1;
    do {
        simpl.emplace_back(in[i]);
        i = parent[i];
    } while (i.valid());
    assert(simpl.back() == in[0]);
    std::reverse(simpl.begin(), simpl.end());
    return Curve("", simpl);
}
