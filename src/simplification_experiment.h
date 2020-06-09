#ifndef SIMPLIFICATION_EXPERIMENT_H
#define SIMPLIFICATION_EXPERIMENT_H

#include <cstddef>
#include <vector>

#include "Curve.h"

namespace experiments {
    std::vector<Curve> sample(std::vector<Curve> const& curves,
        unsigned period);

    void evaluate(std::vector<Curve> const& curves, std::size_t ell);
}
#endif
