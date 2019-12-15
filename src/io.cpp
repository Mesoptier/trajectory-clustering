#include "io.h"

#define FILE_ERROR(filename) throw std::runtime_error("Failed to open file " + filename)

void io::readCurve(const std::string& filename, Curve& curve, int header_size) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        FILE_ERROR(filename);
    }

    distance_t x, y;
    while (file >> x >> y) {
        // Ignore duplicate coordinates
        if (!curve.empty() && approx_equal(curve.back(), {x, y})) {
            continue;
        }
        curve.push_back({x, y});
    }

    file.close();
}

Curve io::readCurve(const std::string& filename, int header_size) {
    Curve curve;
    readCurve(filename, curve, header_size);
    return curve;
}

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