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
     * \return The IDs of new centers.
     */
    std::vector<std::size_t> compute(std::size_t n, std::size_t k,
        DistanceMatrix<distance_t> const& d);

    /**
     * \brief This is an implementation of the PAM algorithm that skips the
     * `build' step and uses the provided initial centers instead.
     * \param n The number of data points.
     * \param d The distance matrix.
     * \param init Initial center IDs. Updated to new center IDs.
     */
    void compute_init(std::size_t n, DistanceMatrix<distance_t> const& d,
        std::vector<std::size_t>& init);
}
#endif
