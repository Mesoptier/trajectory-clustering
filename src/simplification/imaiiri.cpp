#include "simplification/imaiiri.h"

#include <limits>
#include "IntegralFrechet/IntegralFrechet.h"
#include "SymmetricMatrix.h"

namespace {
    struct entry {
        distance_t cost;
        PointID prev;
        entry(): cost(std::numeric_limits<distance_t>::infinity()), prev(0) {}
    };
}

std::pair<distance_t, Curve> simplification::imai_iri::simplify(const Curve& in,
        const PointID& ell) {
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
                distances.at(i, j) = IntegralFrechet(
                    in.slice(i, j), Curve("", {in[i], in[j]}),
                    ParamMetric::L1, 50).compute_matching().cost;
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
            for (PointID j = k; j < i; ++j) {
                distance_t opt = Q[k - 1][j].cost + distances.at(j, i);
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
    return std::make_pair(cost, Curve("", simpl_points));
}
