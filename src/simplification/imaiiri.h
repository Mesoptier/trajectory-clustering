#ifndef IMAIIRI
#define IMAIIRI

#include <utility>
#include "../Curve.h"

namespace simplification::imai_iri {
    /**
     * \brief Run Imai--Iri simplification with CDTW.
     *
     * This does min-epsilon local simplification, producing a polygonal curve
     * with at most ell segments and returns the cost.
     * \param in The original curve.
     * \param ell The maximum complexity of the simplification.
     * \return The cost of the simplification, i.e. the sum of CDTW values over
     * all the simplified segments w.r.t. the matching subcurve, and the
     * simplified curve itself.
     */
    std::pair<distance_t, Curve> simplify(const Curve& in, const PointID& ell);
}
#endif
