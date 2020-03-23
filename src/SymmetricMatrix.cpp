#include "SymmetricMatrix.h"
#include "geom.h"

template<typename T>
SymmetricMatrixT<T> SymmetricMatrixT<T>::read(std::ifstream& file) {
    // Parse header
    std::string line;
    std::getline(file, line);
    if (line != "%%MatrixMarket matrix array real symmetric") {
        throw std::runtime_error("unsupported header");
    }

    // Parse matrix size
    size_t m, n;
    file >> m >> n;
    if (m != n) {
        throw std::runtime_error("matrix is not square");
    }

    // Parse values
    SymmetricMatrixT<T> matrix(n);
    for (auto& value : matrix.data) {
        file >> value;
    }

    return matrix;
}

template<typename T>
void SymmetricMatrixT<T>::write(std::ofstream& file, unsigned int precision) const {
    file << "%%MatrixMarket matrix array real symmetric\n";
    file << n << ' ' << n << '\n';
    file.precision(precision);
    for (T value : data) {
        file << value << '\n';
    }
}

template
class SymmetricMatrixT<distance_t>;
