#pragma once

#include <cstddef>
#include "../SymmetricMatrix.h"

namespace clustering::pam {

    void compute(size_t n, size_t k, const SymmetricMatrix& d);

}