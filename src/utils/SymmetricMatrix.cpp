#include "utils/SymmetricMatrix.h"

#include <stdexcept>
#include <system_error>
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
    std::size_t m, n;
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
void SymmetricMatrixT<T>::write(std::string const& path,
        unsigned int precision) const {
    std::ofstream file(path, std::ios::out | std::ios::trunc);
    if (!file)
        throw std::system_error(errno, std::system_category(),
                                "Failed to open " + path);

    file << "%%MatrixMarket matrix array real symmetric\n";
    file << n << ' ' << n << '\n';
    file.precision(precision);
    for (T value : data) {
        file << value << '\n';
    }

    file.close();
}

template
class SymmetricMatrixT<distance_t>;
