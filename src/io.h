#pragma once

#include "geom.h"
#include "Curve.h"

namespace io {
    void readCurve(const std::string& filename, Curve& curve, int header_size = 0);
    Curve readCurve(const std::string& filename, int header_size = 0);
    void exportPoints(const std::string& filename, const Points& points);
}