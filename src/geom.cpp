#include "geom.h"

distance_t perp(const Point& a, const Point& b) {
    return a(0) * b(1) - a(1) * b(0);
}

bool lessThanMonotone(const Point& a, const Point& b)  {
    return a(0) <= b(0) && a(1) <= b(1);
}

const Point& minMonotone(const Point& a, const Point& b) {
    return lessThanMonotone(a, b) ? a : b;
}

Points makePoints(std::initializer_list<Point> pointsList) {
    Points result(pointsList.size(), 2);
    int row = 0;
    for (const auto& point : pointsList) {
        result.row(row++) = point;
    }
    return result;
}
