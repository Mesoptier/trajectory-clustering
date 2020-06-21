#pragma once

#include <cstddef>
#include "../CurveSimpMatrix.h"


namespace clustering::pam_simp {

    std::vector<size_t> compute(size_t n, size_t k, const CurveSimpMatrix& d, std::vector<size_t>& medoids);
    
}
