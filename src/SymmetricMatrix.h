#ifndef SYMMETRICMATRIXH
#define SYMMETRICMATRIXH

#include <vector>
#include <fstream>

/**
 * Stores a N x N symmetric matrix in a flat array of size (N*(N+1))/2.
 */
template<typename T>
class SymmetricMatrixT
{
public:
    const size_t n;

private:
    std::vector<T> data;

    /**
     * \brief Get the index in the data array from two coordinates.
     * \param i The row (0 to n - 1).
     * \param j The column (0 to n - 1).
     * \return The value in the symmetric matrix at (i, j) or (j, i).
     */
    inline size_t idx(size_t i, size_t j) const {
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
    explicit SymmetricMatrixT(size_t s): n(s), data((n * (n + 1)) / 2, 0) {}

    /**
     * \brief Access the element at (i, j).
     * \param i The row (0 to n - 1).
     * \param j The column (0 to n - 1).
     * \return The (modifiable) element in the matrix.
     */
    T& at(size_t i, size_t j) {
        return data.at(idx(i, j));
    }

    /**
     * \brief Access the element at (i, j) without modification.
     * \param i The row (0 to n - 1).
     * \param j The column (0 to n - 1).
     * \return The element in the matrix.
     */
    const T& at(size_t i, size_t j) const {
        return data.at(idx(i, j));
    }

    /**
     * \brief Write the matrix to a file.
     * \param file The output file.
     * \param precision Floating point precision, if applicable.
     */
    void write(std::ofstream& file, unsigned int precision = 10) const;

    /**
     * \brief Read a matrix from a file.
     * \param file The input file.
     * \return The matrix.
     */
    static SymmetricMatrixT<T> read(std::ifstream& file);
};

extern template class SymmetricMatrixT<double>;
#endif
