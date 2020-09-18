#ifndef REMOVE_STOPS_H
#define REMOVE_STOPS_H

#include "Curve.h"

namespace io {
    /**
     * \brief Remove the stops (small radius parts) in the storks trajectories.
     * \param curve The original curve.
     * \return The curve without stops.
     */
    Curve remove_stops(Curve const& curve);
}
#endif
