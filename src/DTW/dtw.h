#ifndef DTWH
#define DTWH

#include <cstddef>
#include <limits>
#include <utility>
#include "Curve.h"
#include "geom.h"

class DTW {
    // Entry in the dynamic program
    struct entry {
        // cost of the matching so far
        distance_t cost;
        // previous best matching
        std::pair<PointID, PointID> prev;
        /**
         * \brief Initialise the dynamic program entry with infinite cost.
         */
        entry(): cost(std::numeric_limits<distance_t>::infinity()),
                 prev(0, 0) {}
    };
    // The dynamic program
    std::vector<std::vector<entry>> costs;
    // L2 / L1 / Linf for point distance
    Norm const n;
    // The optimal matching of indices
    std::vector<std::pair<PointID, PointID>> m_matching;

    /**
     * \brief Compute the distance between two points with the current metric.
     * \param a The first point.
     * \param b The second point.
     * \return The distance between the two points.
     */
    distance_t dist(Point const& a, Point const& b);

    /**
     * \brief Return the pair of indices that is the best previous step to
     * reach (i, j).
     * \param i The row index in DP; at least 1.
     * \param j The column index in DP; at least 1.
     * \return One of {(i - 1, j - 1), (i - 1, j), (i, j - 1)} that minimises
     * the cost.
     */
    std::pair<std::size_t, std::size_t> min_prev(std::size_t i, std::size_t j);

public:
    /**
     * \brief Compute DTW between two curves using the given distance metric.
     * \param c1 First curve.
     * \param c2 Second curve.
     * \param metric Distance metric.
     */
    DTW(Curve const& c1, Curve const& c2, Norm const metric = Norm::L2);

    /**
     * \brief Return the cost of the optimal matching.
     * \return The DTW cost.
     */
    distance_t cost() const;

    /**
     * \brief Return the optimal matching.
     * \return The optimal matching as a sequence of indices of matched points.
     */
    std::vector<std::pair<PointID, PointID>> const& matching() const;
};
#endif
