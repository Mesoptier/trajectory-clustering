#include "utils/CurveSimpMatrix.h"

#include <fstream>
#include <stdexcept>
#include <system_error>

CurveSimpMatrix::CurveSimpMatrix(std::string const& path) {
    read(path);
}

CurveSimpMatrix::CurveSimpMatrix(std::vector<std::vector<distance_t>> dist):
        matrix(std::move(dist)) {
    // Empty.
}

CurveSimpMatrix::CurveSimpMatrix(Curves const& curves,
        Curves const& simplifications,
        std::function<distance_t(Curve const&, Curve const&)> const& dist)
        : matrix(curves.size(), std::vector<distance_t>(curves.size())) {
    using szt = Curves::size_type;
    if (curves.size() != simplifications.size())
        throw std::runtime_error("There must be the same number of curves and "
            "simplifications in CurveSimpMatrix constructor");

    #pragma omp parallel for collapse(2) schedule(dynamic)
    for (szt i = 0; i < curves.size(); ++i) {
        for (szt j = 0; j < curves.size(); ++j)
            matrix[i][j] = (i == j) ? 0.0 : dist(curves[i], simplifications[j]);
    }
}

void CurveSimpMatrix::write(std::string const& path,
        unsigned int precision) const {
    std::ofstream file(path, std::ios::out | std::ios::trunc);
    if (!file)
        throw std::system_error(errno, std::system_category(),
                                "Failed to open " + path);

    file.precision(precision);
    for (auto const& row: matrix) {
        for (auto const& distance: row) {
            file << distance << " ";
        }
        file << "\n";
    }

    file.close();
}

void CurveSimpMatrix::read(std::string const& path) {
    matrix.clear();

    std::ifstream file(path);
    if (!file)
        throw std::system_error(errno, std::system_category(),
                                "Failed to open " + path);

    std::string line;
    while (std::getline(file, line)) {
        matrix.push_back(std::vector<distance_t>());
        std::istringstream tokenStream(line);
        std::string distance;
        while (std::getline(tokenStream, distance, ' ')) {
            matrix.back().push_back(std::stod(distance));
        }
    }

    file.close();
}

CurveSimpMatrix matrix_from_subset(CurveSimpMatrix const& matrix,
        std::vector<std::size_t> const& indices) {
    std::vector<std::vector<distance_t>> distances(indices.size(),
        std::vector<distance_t>(indices.size()));
    for (unsigned i = 0; i < indices.size(); ++i)
        for (unsigned j = 0; j < indices.size(); ++j)
            distances[i][j] = matrix.at(indices[i], indices[j]);
    return CurveSimpMatrix(distances);
}
