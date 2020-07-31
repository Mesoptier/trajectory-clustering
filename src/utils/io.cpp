#include "utils/io.h"

#include <fstream>
#include <stdexcept>

void io::read_curve(const std::string& filename, Curve& curve,
        std::size_t header_size) {
    std::ifstream file(filename);
    if (!file.is_open())
        throw std::runtime_error("Failed to open file " + filename);

    std::string line;
    while (header_size > 0) {
        std::getline(file, line);
        --header_size;
    }

    distance_t x, y;
    while (file >> x >> y && std::getline(file, line)) {
        // Ignore duplicate coordinates
        if (!curve.empty() && approx_equal(curve.back(), {x, y}))
            continue;
        curve.push_back({x, y});
    }

    file.close();
}

Curve io::read_curve(const std::string& filename, std::size_t header_size) {
    Curve curve(filename);
    read_curve(filename, curve, header_size);
    return curve;
}

std::vector<Curve> io::read_curves(const std::string& directory,
        std::size_t header_size) {
    const auto index_filename = directory + "/dataset.txt";
    std::ifstream index(index_filename);
    if (!index.is_open())
        throw std::runtime_error("Failed to open file " + index_filename);

    std::vector<Curve> curves;

    std::string line;
    while (std::getline(index, line)) {
        auto filename = directory + "/" + line;
        curves.emplace_back(read_curve(filename, header_size));
    }

    return curves;
}

void io::export_points(const std::string& filename, const Points& points) {
    std::ofstream file(filename);
    if (!file.is_open())
        throw std::runtime_error("Failed to open file " + filename);

    file.precision(10);
    for (const auto& point : points)
        file << point.x << ' ' << point.y << '\n';

    file.close();
}
