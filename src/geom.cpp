#include "geom.h"

distance_t perp(const Point& a, const Point& b) {
    return a(0) * b(1) - a(1) * b(0);
}