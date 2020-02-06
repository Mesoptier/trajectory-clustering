#pragma once

#include <vector>
#include <fstream>

/**
 * Stores a N x N symmetric matrix in a flat array of size (N*(N+1))/2.
 */
class SymmetricMatrix
{
public:
    const size_t n;

private:
    std::vector<double> data;

    /**
     * Get the index in the data array from two coordinates.
     *
     * @param i
     * @param j
     * @return
     */
    inline size_t idx(size_t i, size_t j) const {
        if (i <= j) {
            return i * n - (i * (i + 1)) / 2 + j;
        } else {
            return j * n - (j * (j + 1)) / 2 + i;
        }
    }

public:
    explicit SymmetricMatrix(size_t n) : n(n), data((n * (n + 1)) / 2, 0) {}

    /**
     * Access the element at (i, j).
     *
     * @param i
     * @param j
     * @return
     */
    double& at(size_t i, size_t j) {
        return data.at(idx(i, j));
    }
    const double& at(size_t i, size_t j) const {
        return data.at(idx(i, j));
    }

    /**
     * Write the matrix to a file.
     *
     * @param file
     * @param precision
     */
    void write(std::ofstream& file, unsigned int precision = 10) const;

    /**
     * Read a matrix from a file.
     *
     * @param file
     * @return
     */
    static SymmetricMatrix read(std::ifstream& file);

    bool operator==(const SymmetricMatrix& rhs) const {
        return n == rhs.n && data == rhs.data;
    }
};
