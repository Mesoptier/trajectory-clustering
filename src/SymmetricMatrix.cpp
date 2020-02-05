#include "SymmetricMatrix.h"

void SymmetricMatrix::write(std::ofstream& file, unsigned int precision) const {
    file << "%%MatrixMarket matrix array real symmetric\n";
    file << n << ' ' << n << '\n';
    file.precision(precision);
    for (double value : data) {
        file << value << '\n';
    }
}

SymmetricMatrix SymmetricMatrix::read(std::ifstream& file) {
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
    SymmetricMatrix matrix(n);
    for (size_t i = 0; i < n; ++i) {
        file >> matrix.data[i];
    }

    return matrix;
}
