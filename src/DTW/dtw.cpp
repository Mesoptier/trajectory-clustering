#include "dtw.h"

#include <algorithm>
#include <tuple>
#include <vector>

distance_t DTW::dist(const Point& a, const Point& b) {
    auto diff = b - a;
    return norm(diff, n);
}

std::pair<std::size_t, std::size_t> DTW::min_prev(std::size_t i, std::size_t j) {
    distance_t mdist = std::min({costs[i - 1][j].cost, costs[i][j - 1].cost,
                                 costs[i - 1][j - 1].cost});
    if (mdist == costs[i - 1][j - 1].cost)
        return {i - 1, j - 1};
    if (mdist == costs[i - 1][j].cost)
        return {i - 1, j};
    return {i, j - 1};
}

DTW::DTW(const Curve& c1, const Curve& c2, Norm metric):
    costs(c1.size() + 1, std::vector<entry>(c2.size() + 1, entry())),
    n(metric) {
    costs[0][0].cost = 0.0;

    // We use row and column 0 as dummy row / column to avoid special handling.
    for (std::size_t i = 1; i < costs.size(); ++i) {
        for (std::size_t j = 1; j < costs[0].size(); ++j) {
            std::size_t k, l;
            std::tie(k, l) = min_prev(i, j);
            costs[i][j].cost = dist(c1[i - 1], c2[j - 1]) + costs[k][l].cost;
            costs[i][j].prev = {k, l};
        }
    }

    std::size_t i = c1.size(), j = c2.size();
    do {
        m_matching.emplace_back(i - 1, j - 1);
        std::tie(i, j) = costs[i][j].prev;
    } while (i != 0 || j != 0);
    std::reverse(m_matching.begin(), m_matching.end());
}

distance_t DTW::cost() const {
    return costs[costs.size() - 1][costs[0].size() - 1].cost;
}

const std::vector<std::pair<PointID, PointID>>& DTW::matching() const {
    return m_matching;
}
