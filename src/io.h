#pragma once

#include "geom.h"
#include "Curve.h"

namespace io {
    void read_curve(const std::string& filename, Curve& curve, int header_size = 0);
    Curve read_curve(const std::string& filename, int header_size = 0);
    std::vector<Curve> read_curves(const std::string& directory);
    void export_points(const std::string& filename, const Points& points);
    std::vector<Curve> read_pigeon_curves(const std::string& directory);
    std::vector<Curve> read_pigeon_curves_utm(const std::string& directory);
}
