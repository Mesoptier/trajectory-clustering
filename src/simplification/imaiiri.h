#ifndef IMAIIRI
#define IMAIIRI

#include <functional>
#include <utility>
#include "Curve.h"

namespace simplification::imai_iri {
    /**
     * \brief Run Imai--Iri-inspired simplification with a distance function of
     * choice.
     *
     * This does min-epsilon simplification, producing a polygonal curve with at
     * most ell segments, and returns the cost. With max = true, this runs a
     * local simplification (every simplified segment has to be close enough to
     * the subcurve); with max = false, the condition is relaxed and we consider
     * the sum of the costs of individual segments instead. Assuming distance
     * computation takes time proportional to the number of vertices, this runs
     * in time O(ln^3).
     * \param in The original curve.
     * \param ell The maximum complexity of the simplification.
     * \param dist The distance function.
     * \param max Whether to sum of costs over segments or take the largest
     * one-segment cost.
     * \return The cost of the simplification, i.e. the sum / max of dist values
     * over all the simplified segments w.r.t. the matching subcurve, and the
     * simplified curve itself.
     */
    std::pair<distance_t, Curve> simplify(const Curve& in, const PointID& ell,
        const std::function<distance_t(const Curve&, const Curve&)>& dist,
        bool max);

    /**
     * \brief Run Imai--Iri-inspired simplification with a decision solver for a
     * distance function.
     *
     * Note that unlike the method that directly uses dist, here we cannot use
     * the summation, as all we could sum is true / false values. So, this does
     * the traditional Imai--Iri simplification, where for every segment pq and
     * the corresponding subcurve <pq>, it holds that less_than(pq, <pq>, t) =
     * true for some t, and this is the lowest t that yields at most ell
     * segments. Assuming distance decision takes time proportional to the
     * number of vertices, this runs in time O(n^3 log a), where a is the area
     * of the axis-aligned bounding box around the curve.
     * \param in The original curve.
     * \param ell The maximum complexity of the simplification.
     * \param less_than The decision problem solver for a distance function.
     * \return The simplified curve.
     */
    Curve simplify(const Curve& in, const PointID& ell, const
        std::function<bool(const Curve&, const Curve&, distance_t)>& less_than);

    /**
     * \brief Run Imai--Iri-inspired simplification with a decision solver for a
     * distance function (min-ell version for fixed threshold). Assuming
     * distance decision takes time proportional to the number of vertices, this
     * runs in time O(n^3).
     * time 
     * \param in The original curve.
     * \param delta The threshold.
     * \param less_than The decision problem solver for a distance function.
     * \return The simplified curve.
     */
    Curve simplify(const Curve& in, distance_t delta, const
        std::function<bool(const Curve&, const Curve&, distance_t)>& less_than);
}
#endif
