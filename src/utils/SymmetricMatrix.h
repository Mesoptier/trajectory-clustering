#ifndef SYMMETRICMATRIXH
#define SYMMETRICMATRIXH

#include <cstddef>
#include <fstream>
#include <string>
#include <vector>
#include "utils/DistanceMatrix.h"

/**
 * Stores a N x N symmetric matrix in a flat array of size (N*(N+1))/2.
 */
template<typename T>
class SymmetricMatrixT: public DistanceMatrix<T> {
public:
    std::size_t const n;

private:
    std::vector<T> data;

    /**
     * \brief Get the index in the data array from two coordinates.
     * \param i The row (0 to n - 1).
     * \param j The column (0 to n - 1).
     * \return The value in the symmetric matrix at (i, j) or (j, i).
     */
    inline std::size_t idx(std::size_t i, std::size_t j) const {
        if (i <= j) {
            return i * n - (i * (i + 1)) / 2 + j;
        } else {
            return j * n - (j * (j + 1)) / 2 + i;
        }
    }

public:
    /**
     * \brief Initialize symmetric matrix of size s x s.
     * \param s The dimension.
     */
    explicit SymmetricMatrixT(std::size_t s): n(s), data((n * (n + 1)) / 2, 0) {}

    /**
     * \brief Get the number of curves.
     * \return The number of curves.
     */
    std::size_t size() const override {
        return n;
    }

    /**
     * \brief Access the element at (i, j).
     * \param i The row (0 to n - 1).
     * \param j The column (0 to n - 1).
     * \return The (modifiable) element in the matrix.
     */
    T& at(std::size_t i, std::size_t j) override {
        return data.at(idx(i, j));
    }

    /**
     * \brief Access the element at (i, j) without modification.
     * \param i The row (0 to n - 1).
     * \param j The column (0 to n - 1).
     * \return The element in the matrix.
     */
    T const& at(std::size_t i, std::size_t j) const override {
        return data.at(idx(i, j));
    }

    /**
     * \brief Write the matrix to a file.
     * \param path The output file.
     * \param precision Floating point precision, if applicable.
     */
    void write(std::string const& path,
        unsigned int precision = 10) const override;

    /**
     * \brief Read a matrix from a file.
     * \param file The input file.
     * \return The matrix.
     */
    static SymmetricMatrixT<T> read(std::ifstream& file);
};

extern template class SymmetricMatrixT<double>;
#endif
