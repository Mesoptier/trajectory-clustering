#ifndef PAM_H
#define PAM_H

#include <cstddef>
#include <vector>
#include "geom.h"
#include "utils/DistanceMatrix.h"

namespace clustering::pam {
    /**
     * \brief This is an implementation of PAM algorithm (see
     * https://en.wikipedia.org/wiki/K-medoids). This implementation follows
     * https://arxiv.org/pdf/1810.05691.pdf.
     * \param n The number of data points.
     * \param k The number of clusters.
     * \param d The distance / dissimilarity matrix. Does not have to be
     * symmetric: works with both SymmetricMatrix and CurveSimpMatrix.
     */
    std::vector<std::size_t> compute(std::size_t n, std::size_t k,
        DistanceMatrix<distance_t> const& d);
}
#endif
