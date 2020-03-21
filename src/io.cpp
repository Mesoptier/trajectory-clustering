#include <fstream>
#include <iostream>
#include "io.h"

#define FILE_ERROR(filename) throw std::runtime_error("Failed to open file " + filename)

void io::read_curve(const std::string& filename, Curve& curve, int /*header_size */) {
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

Curve io::read_curve(const std::string& filename, int header_size) {
    Curve curve(filename);
    read_curve(filename, curve, header_size);
    return curve;
}

std::vector<Curve> io::read_curves(const std::string& directory) {
    const auto index_filename = directory + "/dataset.txt";
    std::ifstream index(index_filename);

    std::vector<Curve> curves;

    std::string line;
    while (std::getline(index, line)) {
        auto filename = directory + "/" + line;
        curves.push_back(read_curve(filename));
    }

    return curves;
}

Curve read_pigeon_curve(const std::string& filename) {
    // std::cout << filename << "\n";
    Curve curve(filename);
    std::ifstream file(filename);
    std::string line;
    std::getline(file, line);
    std::string x_string, y_string;

    while (std::getline(file, line)) {
        if (line.length() > 10) {
            std::istringstream tokenStream(line);
            std::getline(tokenStream, x_string, '\t');
            std::getline(tokenStream, y_string, '\t');
            distance_t x = std::stod(x_string);
            distance_t y = std::stod(y_string);

            if (!curve.empty() && approx_equal(curve.back(), {x, y})) {
                continue;
            }

            curve.push_back({x, y});
        }
    }
    return curve;
}

std::vector<Curve> io::read_pigeon_curves(const std::string& directory) {
    const auto index_filename = directory + "/dataset.txt";
    std::ifstream index(index_filename);

    std::vector<Curve> curves;

    std::string line;
    while (std::getline(index, line)) {
        auto filename = directory + "/" + line;
        curves.push_back(read_pigeon_curve(filename));
    }

    return curves;
}

void io::export_points(const std::string& filename, const Points& points) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        FILE_ERROR(filename);
    }

    file.precision(10);
    for (const auto& point : points) {
        file << point.x << ' ' << point.y << '\n';
    }

    file.close();
}
