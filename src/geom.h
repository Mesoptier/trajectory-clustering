#ifndef GEOM_H
#define GEOM_H

#include <array>
#include <cassert>
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <vector>

#include "utils/id.h"
#include "utils/SymmetricMatrix.h"
#include "utils/util.h"

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

    Point& operator-=(Point const& point);

    Point operator-(Point const& point) const;

    Point& operator+=(Point const& point);

    Point operator+(Point const& point) const;

    Point& operator*=(distance_t mult);

    Point operator*(distance_t mult) const;

    Point& operator/=(distance_t div);

    Point operator/(distance_t div) const;

    bool operator==(Point const& other) const;

    bool operator!=(Point const& other) const;

    distance_t dist_sqr(Point const& point) const;

    distance_t dist(Point const& point) const;

    friend std::ostream& operator<<(std::ostream& out, Point const& p);
};

template<>
bool approx_equal<Point>(Point const& a, Point const& b, double tol);

/**
 * Computes perp dot product between two vectors.
 * See: http://mathworld.wolfram.com/PerpDotProduct.html
 */
distance_t perp(Point const& a, Point const& b);

distance_t dot(Point const& a, Point const& b);

distance_t acute_angle(Point const& a, Point const& b);

distance_t norm(Point const& point, Norm p = Norm::L2);

Point normalise(Point const& point, Norm p = Norm::L2);

using PointID = ID<Point>;
using Points = std::vector<Point>;

std::ostream& operator<<(std::ostream& out, Points const& points);

//
// Implicit Edge (pair of points)
//
namespace ImplicitEdge {
    /**
     * Get the point that lies dist along the edge.
     */
    Point interpolate_at(Point const& s, Point const& t, distance_t dist);
}

//
// Directions
//
// short for: backward-forward direction
enum class BFDirection {
    Backward = 0,
    Forward = 1,
};

BFDirection getMonotoneDirection(Point const& a, Point const& b);

struct MonotoneComparator {
    BFDirection direction;

    explicit MonotoneComparator(BFDirection dir): direction(dir) {}

    bool operator()(Point const& a, Point const& b);
};

struct Line
{
    Point origin;
    Point direction;

    Line() = default;

    Line(Point const& o, Point const& d)
        : origin(o), direction(normalise(d)) {}

    /**
     * Create line from two points that lie on the line.
     */
    static Line fromTwoPoints(Point const& a, Point const& b) {
        assert(!approx_equal(a, b));
        return {a, b - a};
    }

    /**
     * Create line from a point that lies on the line and its slope.
     */
    static Line fromPointAndSlope(Point p, distance_t slope) {
        if (slope == INFINITY)
            return {p, {0, 1}};
        if (slope == -INFINITY)
            return {p, {0, -1}};
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

    distance_t operator()(Point const& p) const {
        if (std::abs(direction.x) > 0.707)
            return (p.x - origin.x) / (direction.x);
        return (p.y - origin.y) / (direction.y);
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
        if (isVertical())
            return std::numeric_limits<distance_t>::quiet_NaN();
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

    Point closest(Point const& point) const;
    int side(Point const& point) const;

    /**
     * Test whether the given point lies on this line.
     */
    bool includesPoint(Point const& point) const;

    /**
     * Find point at which the two given lines intersect.
     * ASSUMPTION: line1 and line2 are not parallel
     */
    friend Point intersect(Line const& line1, Line const& line2);

    friend bool isParallel(Line const& line1, Line const& line2) {
        return approx_equal(std::abs(dot(line1.direction, line2.direction)),
                            1.0);
    }

    friend bool isSameDirection(Line const& line1, Line const& line2) {
        return approx_equal(dot(line1.direction, line2.direction), 1.0);
    }

    friend bool isOppositeDirection(Line const& line1, Line const& line2) {
        return approx_equal(dot(line1.direction, line2.direction), -1.0);
    }

    /**
     * Test whether the two given lines are perpendicular to each other.
     */
    friend bool isPerpendicular(Line const& line1, Line const& line2) {
        return approx_zero(dot(line1.direction, line2.direction));
    }
};

//
// Data Types for FrechetLight and IntegralFrechet:
//

class CPoint {
private:
    PointID point;
    distance_t fraction;

    void normalize() {
        assert(fraction >= 0.0 && fraction <= 1.0);
        if (approx_equal(fraction, 1.0)) {
            fraction = 0.0;
            ++point;
        }
    }

public:
    CPoint(PointID p, distance_t fr) : point(p), fraction(std::abs(fr)) {
        normalize();
    }

    CPoint() : point(PointID::invalid_value), fraction(0.0) {}

    PointID getPoint() const {
        return point;
    }

    distance_t getFraction() const {
        return fraction;
    }

    distance_t convert() const {
        return static_cast<distance_t>(point) + fraction;
    }

    void setPoint(PointID p) {
        point = p; 
    }

    void setFraction(distance_t frac) {
        fraction = frac;
        normalize();
    }

    bool operator<(CPoint const& other) const {
        return point < other.point ||
            (point == other.point && fraction < other.fraction);
    }

    bool operator<=(CPoint const& other) const {
        return point < other.point ||
            (point == other.point && fraction <= other.fraction);
    }

    bool operator>(CPoint const& other) const {
        return point > other.point ||
            (point == other.point && fraction > other.fraction);
    }

    bool operator>=(CPoint const& other) const {
        return point > other.point ||
            (point == other.point && fraction >= other.fraction);
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
        return point > other || (point == other && fraction > 0.0);
    }

    bool operator<=(PointID other) const {
        return point < other || (point == other && fraction == 0.0);
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
        stream << std::fixed << std::setprecision(10)
            << static_cast<distance_t>(point) + fraction;
        return stream.str();
    }

    friend std::ostream& operator<<(std::ostream& out, CPoint const& p);
};

using CPosition = std::array<CPoint, 2>;
using CPositions = std::vector<CPosition>;

namespace std {
    template<>
    struct hash<CPosition> {
        std::size_t operator()(CPosition const& x) const noexcept {
            std::size_t seed = 0;
            hash_combine(seed, x[0]);
            hash_combine(seed, x[1]);
            return seed;
        }
    };

    template<>
    struct hash<CPoint> {
        std::size_t operator()(CPoint const& x) const noexcept {
            std::size_t seed = 0;
            hash_combine(seed, x.getPoint());
            hash_combine(seed, x.getFraction());
            return seed;
        }
    };
}

std::ostream& operator<<(std::ostream& out, CPosition const& pos);

struct CInterval;

using CIntervals = std::vector<CInterval>;
using CIntervalsID = ID<CIntervals>;
using CIntervalID = std::size_t;
using CIntervalIDs = std::vector<CIntervalID>;

using CPoints = std::vector<CPoint>;

using CPosition = std::array<CPoint, 2>;
using CPositions = std::vector<CPosition>;

using CurveIdx = std::size_t;

//
// Interval
//
struct Interval
{
    distance_t begin;
    distance_t end;

    Interval(): begin(1.), end(0.) {}

    Interval(distance_t b, distance_t e): begin(b), end(e) {}

    bool operator<(Interval const& other) const {
        return begin < other.begin || (begin == other.begin && end < other.end);
    }

    bool is_empty() const {
        return begin > end;
    }

    bool intersects(Interval const& other) const {
        if (is_empty() || other.is_empty())
            return false;
        return (other.begin >= begin && other.begin <= end) ||
            (other.end >= begin && other.end <= end) ||
            (other.begin <= begin && other.end >= end);
    }
};

using Intervals = std::vector<Interval>;

std::ostream& operator<<(std::ostream& out, Interval const& interval);

struct CInterval {
    CPoint begin;
    CPoint end;

    const CInterval* reach_parent = nullptr;
    CPoint fixed {std::numeric_limits<PointID::IDType>::max(), 0.};
    CurveIdx fixed_curve;

    CPosition getLowerRightPos() const { 
        if (fixed_curve == 0)
            return {fixed, begin};
        return {end, fixed};
    }

    CPosition getUpperLeftPos() const { 
        if (fixed_curve == 0)
            return {fixed, end};
        return {begin, fixed};
    }

    CInterval(CPoint b, CPoint e, CPoint f, CurveIdx fc)
        : begin(b), end(e), fixed(f), fixed_curve(fc) {}

    CInterval()
        : begin(std::numeric_limits<PointID::IDType>::max(), 0.),
          end(std::numeric_limits<PointID::IDType>::lowest(), 0.) {}

    CInterval(CPoint const& b, CPoint const& e): begin(b), end(e) {}

    CInterval(PointID point1, distance_t fraction1, PointID point2,
            distance_t fraction2)
        : begin(point1, fraction1), end(point2, fraction2) {}

    CInterval(PointID b, PointID e)
        : begin(b, 0.), end(e, 0.) {}

    bool operator<(CInterval const& other) const {
        return begin < other.begin || (begin == other.begin && end < other.end);
    }

    bool is_empty() const {
        return end < begin;
    }

    void make_empty() {
        begin = {std::numeric_limits<PointID::IDType>::max(), 0.};
        end = {std::numeric_limits<PointID::IDType>::lowest(), 0.};
    }

    void clamp(CPoint const& min, CPoint const& max) {
        begin = std::max(min, begin);
        end = std::min(max, end);
    }
};

std::ostream& operator<<(std::ostream& out, CInterval const& interval);

class IntersectionAlgorithm {
public:
    static constexpr distance_t eps = 1e-8;

   /**
    * Returns which section of the line segment from line_start to line_end is
    * inside the circle given by circle_center and radius. If the circle and
    * line segment do not intersect, the result is the empty Interval (and outer
    * is the empty Interval, too). Otherwise the result is an interval [x,y]
    * such that the distance at x and at y is at most the radius, i.e., [x,y] is
    * a subset of the free interval. The optional output "outer" is an interval
    * strictly containing the free interval. In other words, "outer" is an
    * interval [x',y'] containing [x,y] such that x-x', y'-y <= eps and:
    * If x = 0 then x' = -eps, while if x > 0 then the distance at x' is more
    * than the radius.
    * If y = 1 then y' = 1+eps, while if y < 1 then the distance at y' is more
    * than the radius.
    */
    static Interval intersection_interval(Point circle_center,
        distance_t radius, Point line_start, Point line_end,
        Interval * outer = nullptr);

private:
    IntersectionAlgorithm() {} // Make class static-only

    static inline bool smallDistanceAt(distance_t interpolate, Point line_start,
        Point line_end, Point circle_center, distance_t radius_sqr);

    static inline distance_t distanceAt(distance_t interpolate,
        Point line_start, Point line_end, Point circle_center);

    static constexpr distance_t save_eps = 0.5 * eps;
    static constexpr distance_t save_eps_half = 0.25 * eps;
};

// Ellipse
struct Ellipse {
    Point center;
    distance_t width;
    distance_t height;
    double alpha;

    void invalidate() {
        width = -1.;
        height = -1.;
    }

    bool is_valid() {
        return width >= 0 && height >= 0;
    }
};

Ellipse segmentsToEllipse(Point const& a1, Point const& b1, Point const& a2,
    Point const& b2, distance_t distance);

// Circle
struct Circle {
    Point center;
    distance_t radius;

    Circle() = default;

    Circle(Point const& p, distance_t r) : center(p), radius(r) {}
};

Circle calcMinEnclosingCircle(Points points);

/**
 * \brief Computes the minimum distance from a point to a segment defined by the
 * endpoints source and target.
 * \param source An endpoint of the segment.
 * \param target An endpoint of the segment.
 * \param point The point to measure the distance to.
 * \return The distance from the point to the segment.
 */
distance_t segPointDist(Point const& source, Point const& target,
    Point const& point);
#endif
