#include <iomanip>
#include <iterator>
#include <iostream>
#include "geom.h"

namespace {
    template<typename T>
    T pow2(T d) { return std::pow(d, 2); }
}

/* POINT */

Point& Point::operator-=(const Point& point) {
    x -= point.x;
    y -= point.y;
    return *this;
}

Point Point::operator-(const Point& point) const {
    auto result = *this;
    result -= point;
    return result;
}

Point& Point::operator+=(const Point& point) {
    x += point.x;
    y += point.y;
    return *this;
}

Point Point::operator+(const Point& point) const {
    auto result = *this;
    result += point;
    return result;
}

Point& Point::operator*=(distance_t mult) {
    x *= mult;
    y *= mult;
    return *this;
}

Point Point::operator*(distance_t mult) const {
    auto result = *this;
    result *= mult;
    return result;
}

Point& Point::operator/=(distance_t distance) {
    x /= distance;
    y /= distance;
    return *this;
}

Point Point::operator/(distance_t distance) const {
    return {x / distance, y / distance};
}

bool Point::operator==(const Point& other) const {
    return x == other.x && y == other.y;
}

bool Point::operator!=(const Point& other) const {
    return !(*this == other);
}

distance_t Point::dist_sqr(const Point& point) const {
    return pow2(x - point.x) + pow2(y - point.y);
}

distance_t Point::dist(const Point& point) const {
    return std::sqrt(dist_sqr(point));
}

distance_t norm(const Point& point, Norm p) {
    switch (p) {
        case Norm::L1:
            return std::abs(point.x) + std::abs(point.y);
        case Norm::L2:
            return sqrt(pow2(point.x) + pow2(point.y));
        case Norm::LInf:
            return std::max(std::abs(point.x), std::abs(point.y));
        default:
            throw std::invalid_argument("Unsupported norm");
    }
}

Point normalise(const Point& point, Norm p) {
    return point / norm(point, p);
}

template<>
bool approx_equal<Point>(const Point& a, const Point& b, double tol) {
    return ::approx_equal(a.x, b.x, tol) && ::approx_equal(a.y, b.y, tol);
}

distance_t perp(const Point& a, const Point& b) {
    return a.x * b.y - a.y * b.x;
}

distance_t dot(const Point& a, const Point& b) {
    return a.x * b.x + a.y * b.y;
}

std::ostream& operator<<(std::ostream& out, const Point& p) {
    out << std::setprecision(15)
        << "{" << p.x << ", " << p.y << "}";
    return out;
}

std::ostream& operator<<(std::ostream& out, const Points& points) {
    out << "{";
    auto it = points.begin();
    while (it != points.end()) {
        out << (*it++);
        if (it != points.end()) {
            out << ", ";
        }
    }
    out << "}";
    return out;
}

//
// Implicit Edge (pair of points)
//

namespace ImplicitEdge {

    Point interpolate_at(const Point& s, const Point& t, distance_t dist) {
        const auto len = s.dist(t);

        // Degenerate case
        if (len == 0) {
            return s;
        }

        const auto i = dist / len;
        return s * (1 - i) + t * i;
    }

}

//
// Directions
//

BFDirection getMonotoneDirection(const Point& a, const Point& b) {
#ifndef NDEBUG
    if (approx_equal(a, b)) {
        throw std::logic_error("Points are equal, and therefore not monotone");
    }
#endif

    if (MonotoneComparator(BFDirection::Forward)(a, b)) {
        return BFDirection::Forward;
    }
    if (MonotoneComparator(BFDirection::Backward)(a, b)) {
        return BFDirection::Backward;
    }

    throw std::logic_error("Points are not monotone");
}

bool MonotoneComparator::operator()(const Point& a, const Point& b) {
    return direction == BFDirection::Forward
           ? (a.x < b.x + ABS_TOL && a.y < b.y - ABS_TOL) || (a.x < b.x - ABS_TOL && a.y < b.y + ABS_TOL)
           : (a.x + ABS_TOL > b.x && a.y - ABS_TOL > b.y) || (a.x - ABS_TOL > b.x && a.y + ABS_TOL > b.y);
}

//
// Lines
//

Point Line::closest(const Point& point) const {
    return direction * dot(point - origin, direction) + origin;
}

int Line::side(const Point& point) const {

    if (isVertical()) {
        if (approx_equal(point.x, origin.x))
            return 0;
        if (point.x > origin.x)
            return 1;
        if (point.x < origin.x)
            return -1;
    }

    const auto val = perp(point - origin, direction);
    // std::cout << "val: " << val << std::endl;
    return approx_zero(val) ? 0 : (val > 0 ? 1 : -1);
}

bool Line::includesPoint(const Point& point) const {
    return side(point) == 0;
}

Point intersect(const Line& line1, const Line& line2) {
    const auto t1 = perp(line2.origin - line1.origin, line2.direction) / perp(line1.direction, line2.direction);
    if (!line1.includesPoint(line1(t1)))
        assert(line1.includesPoint(line1(t1)));
    assert(line2.includesPoint(line1(t1)));
    return line1(t1);
}

//
// CPoint
//

std::ostream& operator<<(std::ostream& out, const CPoint& p) {
    out << std::setprecision(15)
        << "(" << static_cast<size_t>(p.point) << " + " << p.fraction << ")";
    return out;
}

std::ostream& operator<<(std::ostream& out, const CPosition& pos) {
    out << "(" << pos[0] << ", " << pos[1] << ")";
    return out;
}

distance_t _l1_dist(Point a, Point b) {
    return fabs(a.x-b.x) + fabs(a.y-b.y); 
}