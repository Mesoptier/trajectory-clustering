#ifndef CURVE_SIMP_MATRIX
#define CURVE_SIMP_MATRIX

#include <cstddef>
#include <functional>
#include <string>
#include <vector>
#include "Curve.h"
#include "utils/DistanceMatrix.h"

using Curves = std::vector<Curve>;

/**
 * \brief Matrix specifically used for storing pairwise distances between curves
 * and simplifications (asymmetric).
 *
 * Suppose we have n curves. Compute  a simplification for each, getting n
 * simplifications. The entry at (i, j) stores the distance between curve i and
 * simplification of curve j. As a special case, the diagonal is all zeros.
 */
class CurveSimpMatrix: public DistanceMatrix<distance_t> {
    std::vector<std::vector<distance_t>> matrix;

public:
    /**
     * \brief Create a new matrix from file (row-major order).
     * \param path The input file.
     */
    explicit CurveSimpMatrix(const std::string& path);

    /**
     * \brief Create a new matrix from data.
     * \param curves The curves.
     * \param simplifications The simplifications, one per curve.
     * \param dist The distance function to use in filling the matrix.
     */
    explicit CurveSimpMatrix(const Curves& curves, const Curves& simplifications,
        const std::function<distance_t(const Curve&, const Curve&)>& dist);

    /**
     * \brief Get the distance between curve i and simplification j.
     * \param i Curve index.
     * \param j Simplification index.
     * \return The distance, or 0 if i = j.
     */
    distance_t& at(std::size_t i, std::size_t j) override {
        return matrix[i][j];
    }

    /**
     * \brief Get the distance between curve i and simplification j.
     * \param i Curve index.
     * \param j Simplification index.
     * \return The distance, or 0 if i = j.
     */
    const distance_t& at(std::size_t i, std::size_t j) const override {
        return matrix[i][j];
    }

    /**
     * \brief Store the matrix in a file, in row-major order.
     * \param path The output file name.
     * \param precision Floating point precision, if applicable.
     */
    void write(const std::string& path,
        unsigned int precision = 10) const override;

    /**
     * \brief Load the matrix from a file, in row-major order, clearing the
     * current contents.
     * \param path The input file name.
     */
    void read(const std::string& path);
};
#endif
