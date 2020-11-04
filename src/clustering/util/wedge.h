#ifndef WEDGE_H
#define WEDGE_H

#include <cstddef>
#include <utility>
#include <vector>
#include "Curve.h"
#include "geom.h"

namespace clustering {
    /**
     * \brief The point matched to a segment of a wedge, with the wedge segment
     * index and the weight of the point.
     */
    struct WedgePoint {
        Point point;
        unsigned matching_segment_index;
        distance_t weight;

        WedgePoint(Point p, unsigned i, distance_t w) : point(std::move(p)),
            matching_segment_index(std::move(i)), weight(std::move(w)) {}
    };
    using WedgePoints = std::vector<WedgePoint>;

    /**
     * \brief The wedge on three consecutive vertices of the center curve with
     * the points of the clustered curves matched to it.
     */
    struct Wedge {
        Points vertices;
        WedgePoints wedge_points;

        Wedge(Points vs, WedgePoints wps) : vertices(std::move(vs)),
            wedge_points(std::move(wps)) {}
    };
    using Wedges = std::vector<Wedge>;

    /**
     * \brief Given a matching, compute the points of curve_2 matched to a
     * segment of curve_1.
     * \param param_space_path The matching.
     * \param curve_1 The center curve, on which the wedge is computed.
     * \param curve_2 One of the curves in the cluster.
     * \param src_index Index of the point in curve_1 marking the start of the
     * segment that we are matching to.
     * \param seg_index 0 or 1, depending on whether the segment starting at
     * src_index is the first or the second segment of the wedge.
     * \return The wedge points matched to the requested segment.
     */
    WedgePoints get_points_matched_to_segment(Points const& param_space_path,
        Curve const& curve_1, Curve const& curve_2,
        std::size_t src_index, unsigned seg_index);

    /**
     * \brief Compute the cost of the wedge (see paper).
     * \param wedge The wedge to evaluate.
     * \return The cost (sum over all the points matched to the wedge).
     */
    distance_t wedge_cost(Wedge const& wedge);

    /**
     * \brief Find the best position for the middle point of the wedge with
     * grid search.
     * \param wedge The wedge, where the middle point should be moved.
     * \param eps The step of the grid.
     * \param radius The number of steps in each direction on the grid.
     * \return The best point position on the grid in terms of wedge cost.
     */
    Point grid_search(Wedge wedge, distance_t eps = 0.125, int radius = 20);
}
#endif
