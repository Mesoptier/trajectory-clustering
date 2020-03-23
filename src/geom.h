#pragma once

#include <cassert>
#include <iomanip>
#include <vector>
#include <sstream>
#include <array>

#include "util.h"
#include "id.h"
#include "SymmetricMatrix.h"

using distance_t = double;
using SymmetricMatrix = SymmetricMatrixT<distance_t>;

enum class Norm {
    L1,
    L2,
    LInf,
};

//
// Point
//

struct Point {
    distance_t x;
    distance_t y;

    Point() = default;
    Point(distance_t dx, distance_t dy) : x(dx), y(dy) {}

    Point& operator-=(const Point& point);
    Point operator-(const Point& point) const;
    Point& operator+=(const Point& point);
    Point operator+(const Point& point) const;
    Point& operator*=(distance_t mult);
    Point operator*(distance_t mult) const;
    Point& operator/=(distance_t distance);
    Point operator/(distance_t distance) const;

    bool operator==(const Point& other) const;
    bool operator!=(const Point& other) const;

    distance_t dist_sqr(const Point& point) const;
    distance_t dist(const Point& point) const;

    friend std::ostream& operator<<(std::ostream& out, const Point& p);
};

template<>
bool approx_equal<Point>(const Point& a, const Point& b, double tol);

/**
 * Computes perp dot product between two vectors.
 * See: http://mathworld.wolfram.com/PerpDotProduct.html
 */
distance_t perp(const Point& a, const Point& b);
distance_t dot(const Point& a, const Point& b);

distance_t norm(const Point& point, Norm p = Norm::L2);
Point normalise(const Point& point, Norm p = Norm::L2);

using PointID = ID<Point>;
using Points = std::vector<Point>;

std::ostream& operator<<(std::ostream& out, const Points& points);

using CurveID = std::size_t;
using CurveIDs = std::vector<CurveID>;

//
// Implicit Edge (pair of points)
//

namespace ImplicitEdge {

    /**
     * Get the point that lies `dist` along the edge.
     */
    Point interpolate_at(const Point& s, const Point& t, distance_t dist);

}

//
// Directions
//

// short for: backward-forward direction
enum class BFDirection {
    Backward = 0,
    Forward = 1,
};
BFDirection getMonotoneDirection(const Point& a, const Point &b);

struct MonotoneComparator {
    BFDirection direction;
    explicit MonotoneComparator(BFDirection dir): direction(dir) {}
    bool operator()(const Point& a, const Point& b);
};

struct Line
{
    Point origin;
    Point direction;

    Line() = default;

    Line(const Point& o, const Point& d)
        : origin(o), direction(normalise(d)) {}

    /**
     * Create line from two points that lie on the line.
     */
    static Line fromTwoPoints(const Point& a, const Point& b) {
        assert(!approx_equal(a, b));
        return {a, b - a};
    }

    /**
     * Create line from a point that lies on the line and its slope.
     */
    static Line fromPointAndSlope(Point p, distance_t slope) {
        if (slope == INFINITY) {
            return {p, {0, 1}};
        }
        if (slope == -INFINITY) {
            return {p, {0, -1}};
        }
        return {p, {1, slope}};
    }

    static Line horizontal(Point p) {
        return {p, {1, 0}};
    }

    static Line vertical(Point p) {
        return {p, {0, 1}};
    }

    Point operator()(distance_t t) const {
        return origin + direction * t;
    }

    distance_t operator()(const Point& p) const {
        if (std::abs(direction.x) > 0.707) {
            return (p.x - origin.x) / (direction.x);
        } else {
            return (p.y - origin.y) / (direction.y);
        }
    }

    inline bool isVertical() const {
        return approx_zero(direction.x);
    }

    inline bool isHorizontal() const {
        return approx_zero(direction.y);
    }

    /**
     * Get Y coordinate for given X coordinate.
     */
    inline distance_t getY(distance_t x) const {
        if (isVertical()) {
            return NAN;
        }
        return (x - origin.x) * direction.y / direction.x + origin.y;
    }

    /**
     * Get X coordinate for given Y coordinate.
     */
    inline distance_t getX(distance_t y) const {
        if (isHorizontal()) {
            return NAN;
        }
        return (y - origin.y) * direction.x / direction.y + origin.x;
    }

    Point closest(const Point& point) const;
    int side(const Point& point) const;

    /**
     * Test whether the given point lies on this line.
     */
    bool includesPoint(const Point& point) const;

    /**
     * Find point at which the two given lines intersect.
     * ASSUMPTION: line1 and line2 are not parallel
     */
    friend Point intersect(const Line& line1, const Line& line2);

    friend bool isParallel(const Line& line1, const Line& line2) {
        return approx_equal(std::abs(dot(line1.direction, line2.direction)), 1.0);
    }

    friend bool isSameDirection(const Line& line1, const Line& line2) {
        return approx_equal(dot(line1.direction, line2.direction), 1.0);
    }

    friend bool isOppositeDirection(const Line& line1, const Line& line2) {
        return approx_equal(dot(line1.direction, line2.direction), -1.0);
    }

    /**
     * Test whether the two given lines are perpendicular to each other.
     */
    friend bool isPerpendicular(const Line& line1, const Line& line2) {
        return approx_zero(dot(line1.direction, line2.direction));
    }
};


//
// Data Types for FrechetLight and IntegralFrechet:
//

class CPoint
{
private:
    PointID point;
    distance_t fraction;

    void normalize() {
        assert(fraction >= 0. && fraction <= 1.);
        if (fraction == 1.) {
            fraction = 0.;
            ++point;
        }
    }
public:
    CPoint(PointID p, distance_t fr) : point(p), fraction(fr) {
        normalize();
    }
    CPoint() : point(PointID::invalid_value), fraction(0.) {}
    // CPoint(CPoint const&) = default;

    PointID getPoint() const { return point; }
    distance_t getFraction() const { return fraction; }
    distance_t convert() const { return static_cast<distance_t>(point) + fraction; }
    void setPoint(PointID p) { point = p; }
    void setFraction(distance_t frac) {
        fraction = frac;
        normalize();
    }

    bool operator<(CPoint const& other) const {
        return point < other.point || (point == other.point && fraction < other.fraction);
    }
    bool operator<=(CPoint const& other) const {
        return point < other.point || (point == other.point && fraction <= other.fraction);
    }
    bool operator>(CPoint const& other) const {
        return point > other.point || (point == other.point && fraction > other.fraction);
    }
    bool operator>=(CPoint const& other) const {
        return point > other.point || (point == other.point && fraction >= other.fraction);
    }
    bool operator==(CPoint const& other) const {
        return point == other.point && fraction == other.fraction;
    }
    bool operator!=(CPoint const& other) const {
        return point != other.point || fraction != other.fraction;
    }
    bool operator<(PointID other) const {
        return point < other;
    }
    bool operator>(PointID other) const {
        return point > other || (point == other && fraction > 0.);
    }
    bool operator<=(PointID other) const {
        return point < other || (point == other && fraction == 0.);
    }
    bool operator>=(PointID other) const {
        return point >= other;
    }
    bool operator==(PointID other) const {
        return point == other && fraction == 0.;
    }
    bool operator!=(PointID other) const {
        return !(point == other);
    }
    CPoint operator+(distance_t other) const {
        assert(other <= 1.);
        PointID p = point;
        distance_t f = fraction + other;
        if (f > 1.) {
            ++p;
            f -= 1.;
        }
        return CPoint(p, f);
    }
    CPoint operator-(distance_t other) const {
        assert(other <= 1.);
        PointID p = point;
        distance_t f = fraction - other;
        if (f < 0.) {
            --p;
            f += 1.;
        }
        return CPoint(p, f);
    }
    CPoint ceil() const {
        return fraction > 0 ? CPoint(point + 1, 0.) : CPoint(point, 0.);
    }
    CPoint floor() const {
        return CPoint(point, 0.);
    }
    std::string to_string() const {
        std::stringstream stream;
        stream << std::fixed << std::setprecision(10) << static_cast<distance_t>(point) + fraction;
        return stream.str();
    }

    friend std::ostream& operator<<(std::ostream& out, const CPoint& p);
};

using CPosition = std::array<CPoint, 2>;
using CPositions = std::vector<CPosition>;

namespace std {
    template<>
    struct hash<CPosition> {
        size_t operator()(const CPosition& x) const {
            size_t seed = 0;
            hash_combine(seed, x[0]);
            hash_combine(seed, x[1]);
            return seed;
        }
    };

    template<>
    struct hash<CPoint> {
        size_t operator()(const CPoint& x) const {
            size_t seed = 0;
            hash_combine(seed, x.getPoint());
            hash_combine(seed, x.getFraction());
            return seed;
        }
    };
}

std::ostream& operator<<(std::ostream& out, const CPosition& pos);
