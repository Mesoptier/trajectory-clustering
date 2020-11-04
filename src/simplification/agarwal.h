#ifndef AGARWAL_GREEDY
#define AGARWAL_GREEDY

#include <functional>
#include "Curve.h"

namespace simplification::greedy {
    /**
     * \brief Run local simplification approximation by Agarwal et al. with the
     * decision version of a distance function (fixed l, min delta).
     * \param in The original curve.
     * \param ell The maximum complexity of the simplification.
     * \param less_than The decision version of the distance function between
     * two curves, returns true if d(c1, c2) < dist.
     * \return The simplified curve.
     */
    Curve simplify(Curve const& in, std::size_t ell,
        std::function<bool(Curve const&, Curve const&, distance_t)> const&
        less_than);

    /**
     * \brief Run local simplification approximation by Agarwal et al. with the
     * decision version of a distance function (fixed delta, min l).
     * \param in The original curve.
     * \param delta The distance threshold.
     * \param less_than The decision version of the distance function between
     * two curves, returns true if d(c1, c2) < dist.
     * \return The simplified curve.
     */
    Curve simplify(Curve const& in, distance_t delta,
        std::function<bool(Curve const&, Curve const&, distance_t)> const&
        less_than);
}
#endif
