#include "utils/CurveSimpMatrix.h"

#include <fstream>
#include <stdexcept>
#include <system_error>

CurveSimpMatrix::CurveSimpMatrix(const std::string& path) {
    read(path);
}

CurveSimpMatrix::CurveSimpMatrix(const Curves& curves, const Curves&
        simplifications, const std::function<distance_t(const Curve&,
        const Curve&)>& dist) : matrix(curves.size()) {
    using szt = Curves::size_type;
    if (curves.size() != simplifications.size())
        throw std::runtime_error("There must be the same number of curves and "
            "simplifications in CurveSimpMatrix constructor");

    for (szt i = 0; i < curves.size(); ++i) {
        for (szt j = 0; j < curves.size(); ++j) {
            if (i == j)
                matrix[i].push_back(0);
            else
                matrix[i].push_back(dist(curves[i], simplifications[j]));
        }
    }
}

void CurveSimpMatrix::write(const std::string& path,
        unsigned int precision) const {
    std::ofstream file(path, std::ios::out | std::ios::trunc);
    if (!file)
        throw std::system_error(errno, std::system_category(),
                                "Failed to open " + path);

    file.precision(precision);
    for (const auto& row: matrix) {
        for (const auto& distance: row) {
            file << distance << " ";
        }
        file << "\n";
    }

    file.close();
}

void CurveSimpMatrix::read(const std::string& path) {
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
