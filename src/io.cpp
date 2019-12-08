#include "io.h"

#define FILE_ERROR(filename) throw std::runtime_error("Failed to open file " + filename)

void io::exportPoints(const std::string& filename, const Points& points) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        FILE_ERROR(filename);
    }

    file.precision(10);
    for (const auto& point : points) {
        file << point.x << ',' << point.y << '\n';
    }

    file.close();
}
