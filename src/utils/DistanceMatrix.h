#ifndef DISTANCE_MATRIX_H
#define DISTANCE_MATRIX_H

#include <cstddef>
#include <string>

template<typename T>
class DistanceMatrix {
public:
    // The following is to avoid -Wdeprecated warnings of having a destructor
    // but not copy / move constructor / assignment.
    DistanceMatrix() = default;
    DistanceMatrix(DistanceMatrix<T> const&) = default;
    DistanceMatrix(DistanceMatrix<T>&&) noexcept = default;
    DistanceMatrix<T>& operator=(DistanceMatrix<T> const&) = default;
    DistanceMatrix<T>& operator=(DistanceMatrix<T>&&) noexcept = default;
    virtual ~DistanceMatrix() = default;

    /**
     * \brief Access the element at (i, j).
     * \param i The row (0 to n - 1).
     * \param j The column (0 to n - 1).
     * \return The (modifiable) element in the matrix.
     */
    virtual T& at(std::size_t i, std::size_t j) = 0;

    /**
     * \brief Access the element at (i, j) without modification.
     * \param i The row (0 to n - 1).
     * \param j The column (0 to n - 1).
     * \return The element in the matrix.
     */
    virtual const T& at(std::size_t i, std::size_t j) const = 0;

    /**
     * \brief Write the matrix to a file.
     * \param path The output file.
     * \param precision Floating point precision, if applicable.
     */
    virtual void write(const std::string& path,
        unsigned int precision = 10) const = 0;
};
#endif
