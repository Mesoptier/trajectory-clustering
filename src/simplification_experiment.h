#ifndef SIMPLIFICATION_EXPERIMENT_H
#define SIMPLIFICATION_EXPERIMENT_H

#include <cstddef>
#include <vector>

#include "Curve.h"

/**
 * \brief The experiments set-up.
 */
namespace experiments {
    /**
     * \brief Sample a curve every period entries and subsample it to length
     * points.
     * \param curves The full dataset.
     * \param period Select curve i iff i % period == 0.
     * \param length Take first, last, and every kth point of each curve, where
     * k = curve.size() / length.
     * \return The sampled curves.
     */
    std::vector<Curve> sample(std::vector<Curve> const& curves,
        unsigned period, unsigned length);

    /**
     * \brief Run the experiment: simplify each curve to complexity ell with
     * Imai--Iri and greedy approaches using Frechet, (fast) CDTW, and DTW.
     * Report results in a table. Save the worst simplifications and curves to
     * files and show some statistics.
     * \param curves The curves to simplify.
     * \param ell The final complexity of simplifications.
     * \param prefix The files are saved in ./simpl/prefix.
     * \param scale Estimated number of digits before decimal separator for the
     * distance values. Used for pretty printing.
     */
    void evaluate(std::vector<Curve> const& curves, std::size_t ell,
        std::string const& prefix, unsigned scale = 4);
}
#endif
