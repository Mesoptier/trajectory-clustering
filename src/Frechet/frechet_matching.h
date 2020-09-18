#ifndef FRECHET_MATCHING_H
#define FRECHET_MATCHING_H

#include "Curve.h"
#include "geom.h"

namespace frechet {
    Points calcMatching(Curve const& curve1, Curve const& curve2);
}
#endif
