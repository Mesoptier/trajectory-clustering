#include <iomanip>
#include <iterator>
#include <iostream>
#include "geom.h"
#define PI 3.14159265

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

distance_t acute_angle(const Point& a, const Point& b) {
	return std::acos(
		dot(a, b) / (norm(a) * norm(b))
	) / PI;
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
    const auto val = perp(point - origin, direction);
    return approx_zero(val) ? 0 : (val > 0 ? 1 : -1);
}

bool Line::includesPoint(const Point& point) const {
    return side(point) == 0;
}

Point intersect(const Line& line1, const Line& line2) {
	// std::cout << line1.origin << "\n";
	// std::cout << line1.direction << "\n";
	// std::cout << line2.origin << "\n";
	// std::cout << line2.direction << "\n";
    const auto t1 = perp(line2.origin - line1.origin, line2.direction) / perp(line1.direction, line2.direction);
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


//
// intersection_interval
//

inline bool IntersectionAlgorithm::smallDistanceAt(distance_t interpolate, Point line_start, Point line_end, Point circle_center, distance_t radius_sqr) {
	return circle_center.dist_sqr(line_start * (1. - interpolate) + line_end * interpolate) <= radius_sqr;
}

inline distance_t IntersectionAlgorithm::distanceAt(distance_t interpolate, Point line_start, Point line_end, Point circle_center) {
	return circle_center.dist_sqr(line_start * (1. - interpolate) + line_end * interpolate);
}

Interval IntersectionAlgorithm::intersection_interval(Point circle_center, distance_t radius, Point line_start, Point line_end, Interval * outer /* = nullptr*/)
{
    // The line can be represented as line_start + lambda * v
    const Point v = line_end - line_start;
	const distance_t rad_sqr = radius * radius;
	
    // Find points p = line_start + lambda * v with
    //     dist(p, circle_center) = radius
    // <=> sqrt(p.x^2 + p.y^2) = radius
    // <=> p.x^2 + p.y^2 = radius^2
    // <=> (line_start.x + lambda * v.x)^2 + (line_start.y + lambda * v.y)^2 = radius^2
    // <=> (line_start.x^2 + 2 * line_start.x * lambda * v.x + lambda^2 * v.x^2) + (line_start.y^2 + 2 * line_start.y * lambda * v.y + lambda^2 * v.y^2) = radius^2
    // <=> lambda^2 * (v.x^2 + v.y^2) + lambda * (2 * line_start.x * v.x + 2 * line_start.y * v.y) + line_start.x^2 + line_start.y^2) - radius^2 = 0
    // let a := v.x^2 + v.y^2, 
	// let b := line_start.x * v.x + line_start.y * v.y, 
	// let c := line_start.x^2 + line_start.y^2 - radius^2
    // <=> lambda^2 * a + lambda * 2 b + c = 0
    // <=> lambda^2 + (2 b / a) * lambda + (c / a) = 0
    // <=> lambda1/2 = - (b / a) +/- sqrt((b / a)^2 - c / a)
	
    const distance_t a = pow2(v.x) + pow2(v.y);
    const distance_t b = (line_start.x - circle_center.x) * v.x + (line_start.y - circle_center.y) * v.y;
    const distance_t c = pow2(line_start.x - circle_center.x) + pow2(line_start.y - circle_center.y) - pow2(radius);

	distance_t mid = - b / a;
    distance_t discriminant = pow2(mid) - c / a;

	const bool smallDistAtZero = smallDistanceAt(0., line_start, line_end, circle_center, rad_sqr);
	const bool smallDistAtOne = smallDistanceAt(1., line_start, line_end, circle_center, rad_sqr);
	bool smallDistAtMid = smallDistanceAt(mid, line_start, line_end, circle_center, rad_sqr);
	
	if (smallDistAtZero && smallDistAtOne) {
		if (outer != nullptr) { *outer = Interval(-eps, 1. + eps); }
		return Interval(0, 1);
	}
	
	if (!smallDistAtMid && smallDistAtZero) {
		mid = 0.;
		smallDistAtMid = true;
	}
	else if (!smallDistAtMid && smallDistAtOne) {
		mid = 1.;
		smallDistAtMid = true;
	}
	
	// Here we need the guarantee that if the free interval has length at least eps
	// then at mid the distance is <=radius
	// This is an assumption about the precision of distance_t computations
	// All remaining rules are free of such assumptions! 
	// (except for trivial ones like this: x + y and x - y have distance at most 2y up to negligible error)
    if (!smallDistAtMid) {
		if (outer != nullptr) { *outer = Interval(); }
		return Interval(); // no intersection;
    }
	
	if (mid <= 0. and !smallDistAtZero) {
		if (outer != nullptr) { *outer = Interval(); }
		return Interval();
	}
	if (mid >= 1. and !smallDistAtOne) {
		if (outer != nullptr) { *outer = Interval(); }
		return Interval();
	}
	
	discriminant = std::max<distance_t>(discriminant, 0.);
	distance_t sqrt_discr = 0.;
	bool sqrt_discr_computed = false;
	distance_t begin, end;
	
	if (smallDistAtZero) {
		begin = 0.;
		if (outer != nullptr) { outer->begin = -eps; }
	}
	else {
		sqrt_discr = std::sqrt(discriminant);
		sqrt_discr_computed = true;
		
		const distance_t lambda1 = mid - sqrt_discr;
		const distance_t innershift = std::min<distance_t>(lambda1 + save_eps_half, std::min<distance_t>(1., mid));
		const distance_t outershift = lambda1 - save_eps_half;
		if (innershift >= outershift && smallDistanceAt(innershift, line_start, line_end, circle_center, rad_sqr) && !smallDistanceAt(outershift, line_start, line_end, circle_center, rad_sqr)) {
			begin = innershift;
			if (outer != nullptr) { outer->begin = outershift; }
		}
		else {
			distance_t left = 0., right = std::min<distance_t>(mid, 1.);
			// invariants throughout binary search:
			//  * !smallDistanceAt(left)
			//  * smallDistanceAt(right)
			//  * 0 <= left <= right <= min(mid,1)
			// Clearly this is stays true after an iteration.
			// Why is it true in the beginning?
			// If smallDistanceAt(0.) then begin would already be set (fourth rule).
			// If !smallDistanceAt(right), then either !smallDistanceAt(mid), contradicting the very first rule, 
			//  or mid >= 1. and smallDistanceAt(1.), contradicting the third rule.
			// Finally, since !smallDistanceAt(left) we cannot have mid <= 0 by the second rule. Thus, right = min(mid,1) >= 0. = left
			while (right - left > save_eps) {
				distance_t m = 0.5 * (left + right);
				if (smallDistanceAt(m, line_start, line_end, circle_center, rad_sqr)) { right = m; }
				else { left = m; }
			}
			begin = right;
			if (outer != nullptr) { outer->begin = left; }
		}
	}
	
	if (smallDistAtOne) {
		end = 1.;
		if (outer != nullptr) { outer->end = 1. + eps; }
	}
	else {
		if (!sqrt_discr_computed) {
			sqrt_discr = std::sqrt(discriminant);
		}
		
		const distance_t lambda2 = mid + sqrt_discr;
		const distance_t innershift = std::max<distance_t>(lambda2 - save_eps_half, std::max<distance_t>(0., mid));
		const distance_t outershift = lambda2 + save_eps_half;
		if (innershift <= outershift && smallDistanceAt(innershift, line_start, line_end, circle_center, rad_sqr) && !smallDistanceAt(outershift, line_start, line_end, circle_center, rad_sqr)) {
			end = innershift;
			if (outer != nullptr) { outer->end = outershift; }
		}
		else {
			distance_t left = std::max<distance_t>(mid, 0.), right = 1.;
			// invariants throughout binary search:
			//  * smallDistanceAt(left)
			//  * !smallDistanceAt(right)
			//  * max(mid,0) <= left <= right <= 1
			while (right - left > save_eps) {
				distance_t m = 0.5 * (left + right);
				if (smallDistanceAt(m, line_start, line_end, circle_center, rad_sqr)) { left = m; }
				else { right = m; }
			}
			end = left;
			if (outer != nullptr) { outer->end = right; }
		}
	}
	
	assert(smallDistanceAt(begin, line_start, line_end, circle_center, rad_sqr));
	assert(smallDistanceAt(end, line_start, line_end, circle_center, rad_sqr));
	assert(0. <= begin && begin <= end && end <= 1.);
	
	assert(outer == nullptr || outer->begin < 0. || !smallDistanceAt(outer->begin, line_start, line_end, circle_center, rad_sqr));
	assert(outer == nullptr || outer->end > 1. || !smallDistanceAt(outer->end, line_start, line_end, circle_center, rad_sqr));
	assert(outer == nullptr || (outer->begin < begin && begin - outer->begin <= eps));
	assert(outer == nullptr || (outer->end > end && outer->end - end <= eps));
	
	assert(outer == nullptr || (outer->begin <= 1.));
	assert(outer == nullptr || (outer->end >= 0.));
	
	// These asssertions fail - use exact arithmetic to make it work??
	//assert(begin - eps < 0. || !smallDistanceAt(begin - eps, line_start, line_end, circle_center, rad_sqr));
	//assert(end + eps > 1. || !smallDistanceAt(end + eps, line_start, line_end, circle_center, rad_sqr));
	
	return Interval{ begin, end };
}

Ellipse segmentsToEllipse(Point const& a1, Point const& b1, Point const& a2, Point const& b2, distance_t distance)
{
	Ellipse e;

	auto const pi = std::atan2(0, -1);

	// Check if segments are parallel
	auto dir1 = b1 - a1;
	auto dir2 = b2 - a2;
	dir1 /= (std::sqrt(pow2(dir1.x)+pow2(dir1.y)));
	dir2 /= (std::sqrt(pow2(dir2.x)+pow2(dir2.y)));
	if (std::abs(dir1.x*dir2.x + dir1.y*dir2.y) >= 0.999) {
		e.invalidate();
		return e;
	}

	// First assign coefficients of x^2, y^2, ...
	// Those are calculated by hand
	auto A = pow2(a1.x) - 2*a1.x*b1.x + pow2(a1.y) - 2*a1.y*b1.y + pow2(b1.x) + pow2(b1.y);
	auto B = -2*a1.x*a2.x + 2*a1.x*b2.x - 2*a1.y*a2.y + 2*a1.y*b2.y + 2*a2.x*b1.x + 2*a2.y*b1.y - 2*b1.x*b2.x - 2*b1.y*b2.y;
	auto C = pow2(a2.x) - 2*a2.x*b2.x + pow2(a2.y) - 2*a2.y*b2.y + pow2(b2.x) + pow2(b2.y);
	auto D = 2*a1.x*b1.x - 2*a1.x*b2.x + 2*a1.y*b1.y - 2*a1.y*b2.y - 2*pow2(b1.x) + 2*b1.x*b2.x - 2*pow2(b1.y) + 2*b1.y*b2.y;
	auto E = -2*a2.x*b1.x + 2*a2.x*b2.x - 2*a2.y*b1.y + 2*a2.y*b2.y + 2*b1.x*b2.x + 2*b1.y*b2.y - 2*pow2(b2.x) - 2*pow2(b2.y);
	auto F = pow2(b1.x) - 2*b1.x*b2.x + pow2(b1.y) - 2*b1.y*b2.y + pow2(b2.x) + pow2(b2.y) - pow2(distance);

	// This should not fail if they are not parallel
	assert(pow2(B) - 4*A*C <= 0.0);

	// Now convert those to the form of the Ellipse type
	// See: https://en.wikipedia.org/wiki/Ellipse#Canonical_form
	e.center.x = (2*C*D - B*E)/(pow2(B) - 4*A*C);
	e.center.y = (2*A*E - B*D)/(pow2(B) - 4*A*C);
	e.width = -std::sqrt(2*(A*pow2(E) + C*pow2(D) - B*D*E + (pow2(B) - 4*A*C)*F)*(A+C+std::sqrt(pow2(A-C) + pow2(B))))/(pow2(B) - 4*A*C);
	e.height = -std::sqrt(2*(A*pow2(E) + C*pow2(D) - B*D*E + (pow2(B) - 4*A*C)*F)*(A+C-std::sqrt(pow2(A-C) + pow2(B))))/(pow2(B) - 4*A*C);

	if (B == 0) {
		e.alpha = A < C ? 0 : pi/2;
	}
	else {
		e.alpha = std::atan((C-A-std::sqrt(pow2(A-C) + pow2(B)))/B);
	}
	e.alpha = e.alpha*360/(2*pi);

	return e;
}


namespace
{

bool inCircle(Point const& p, Circle const& circle) {
	const distance_t eps = 0.00000000000001;
	return p.dist_sqr(circle.center) <= circle.radius + eps;
}

Circle circleFromTwo(Point p1, Point p2)
{
	auto center  = (p1 + p2)/2.;
	auto radius = p1.dist_sqr(center);
	return Circle(center, radius);
}

Circle calcCircle(Point p1, Point p2, Point p3)
{
	static constexpr distance_t epsilon = 1e-10;

	// compute normal vectors and their offset vectors
	auto q1 = (p1 + p2)/2.;
	auto q2 = (p2 + p3)/2.;
	auto n1 = Point((p2 - p1).y, -1*(p2 - p1).x);
	auto n2 = Point((p3 - p2).y, -1*(p3 - p2).x);

	// Avoid division by zero. Due to non-colinearity n1.x is non-zero after the swap.
	if (n1.x == 0) {
		std::swap(n1, n2);
		std::swap(q1, q2);
	}
	assert(n1.x != 0);

	// main calculation
	auto a = q1.y - q2.y + (q2.x - q1.x)*n1.y/n1.x;
	auto b = n2.y - (n2.x*n1.y)/n1.x;

	// if the points p1, p2, p3 are colinear
	if (std::abs(b) < epsilon) {
		// rather naive, but well...
		if (inCircle(p1, circleFromTwo(p2,p3))) { return circleFromTwo(p2,p3); }
		if (inCircle(p2, circleFromTwo(p1,p3))) { return circleFromTwo(p1,p3); }
		return circleFromTwo(p1,p2);
	}

	auto i2 = a/b;

	auto center = q2 + n2*i2;
	auto radius = center.dist_sqr(p1);
	return Circle(center, radius);
}

} // end anonymous namespace


Circle calcMinEnclosingCircle(Points points)
{
	if (points.size() == 0) {
		return Circle();
	}
	if (points.size() == 1) {
		return Circle(points.front(), 0.);
	}

	std::random_shuffle(points.begin(), points.end());

	auto current_circle = circleFromTwo(points[0], points[1]);
	for (PointID i = 2; i < points.size(); ++i) {
		auto p_i = points[i];
		if (inCircle(p_i, current_circle)) { continue; }

		// from here: p_i not in current_circle
		current_circle = circleFromTwo(p_i, points[0]);
		for (PointID j = 1; j < i; ++j) {
			auto p_j = points[j];
			if (inCircle(p_j, current_circle)) { continue; }

			// from here: p_j not in current_circle
			current_circle = circleFromTwo(p_i, p_j);
			for (PointID k = 0; k < j; ++k) {
				auto p_k = points[k];
				if (inCircle(p_k, current_circle)) { continue; }

				// from here: p_k not in current_circle
				current_circle = calcCircle(p_i, p_j, p_k);
			}
		}
	}

	current_circle.radius = std::sqrt(current_circle.radius);
	return current_circle;
}

// Computes the minimum distance from a point to a segment
// defined by endpoints source and target
distance_t segPointDist(Point& source, Point& target, Point& point) {
	Line line = Line::fromTwoPoints(source, target);
	Point line_vec = target - source;
	Point normal_vec = {-line_vec.y, line_vec.x};

	Line normal_src = Line(source, normal_vec);
	Line normal_tgt = Line(target, normal_vec);

	if (normal_src.side(point) != normal_tgt.side(point)) 
		return point.dist(line.closest(point));
	
	return std::min(point.dist(source), point.dist(target));
}