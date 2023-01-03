#include "2d_utils.h"


namespace utils {


long double angle(Point s, Point t) {
    double PI = 3.14159265359;

    if (approx_equal(s.x, t.x)) {
        if (s.y > t.y)
            return 3*PI/2;
        else
            return PI/2;
    }

    long double acute_angle = atan(fabs(s.y - t.y)/fabs(s.x - t.x));

    if (t.x >= s.x) {
        if (t.y >= s.y)
            return acute_angle;
        else
            return 2*PI - acute_angle;
    } else {
        if (t.y >= s.y)
            return PI - acute_angle;
        else
            return PI + acute_angle;
    }
}

double _sin(double angle) {
    double PI = 3.14159265359;
    if (approx_zero(angle))
        return 0;

    if (approx_equal(angle, PI)) {
        return 0;
    }

    if (approx_equal(angle, 2*PI))
        return 0;

    return sin(angle);
}

double _cos(double angle) {
    double PI = 3.14159265359;
    if (approx_equal(angle, PI/2))
        return 0;

    if (approx_equal(angle, 3*PI/2)) {
        return 0;
    }

    return cos(angle);
}

bool domain_covered(std::vector<ConstrainedBivariatePolynomial<2>> costs, const Cell& cell) {
    std::vector<Constraint> constraints = std::vector<Constraint>();

    double sx = cell.s.x;
    double sy = cell.s.y;
    double tx = cell.t.x;
    double ty = cell.t.y;

    for (auto cost: costs) {
        for (auto constr: cost.left_constraints) {
            constraints.push_back(
                Constraint(false, true, 0, constr)
            );
        }

        for (auto constr: cost.right_constraints) {
            constraints.push_back(
                Constraint(false, false, 0, constr)
            );
        }

        constraints.push_back(Constraint(true, true, cost.y_interval.min, Polynomial<1>()));
        constraints.push_back(Constraint(true, false, cost.y_interval.max, Polynomial<1>()));
    }

    // Check to see if each corner of the domin appears in the
    // feasible region of a cost polynomial

    bool bl = false;
    bool br = false;
    bool tl = false;
    bool tr = false;

    for (auto constr: constraints) {
        bl = bl || constr.valid_point({sx, sy});
        br = br || constr.valid_point({tx, sy});
        tl = tl || constr.valid_point({sx, ty});
        tr = tr || constr.valid_point({tx, ty});
    }

    // if (!bl)
    //     std::cout << "bottom left not covered...\n";
    // if (!br)
    //     std::cout << "bottom right not covered...\n";
    // if (!tl)
    //     std::cout << "top left not covered...\n";
    // if (!tr)
    //     std::cout << "top right not covered...\n";


    return bl && br && tl && tr;
}

/**
 * @brief checks whether a and b are parallel
 * 
 * @param a 
 * @param b 
 * @return true 
 * @return false 
 */
bool is_parallel(Constraint a, Constraint b) {
    if (a.is_y_constraint && b.is_y_constraint)
        return true;

    if (!a.is_y_constraint && !b.is_y_constraint) {
        return approx_equal(a.x_bound.coefficients[1], b.x_bound.coefficients[1]);
    }

    return false;
}

Point intersect(Constraint a, Constraint b) {
    assert(!is_parallel(a, b));

    if (a.is_y_constraint) {
        assert(!b.is_y_constraint);
        double x = b.x_bound.coefficients[1] * a.y_bound + b.x_bound.coefficients[0];
        return {x, a.y_bound};
    }

    if (b.is_y_constraint) {
        assert(!a.is_y_constraint);
        double x = a.x_bound.coefficients[1] * b.y_bound + a.x_bound.coefficients[0];
        return {x, b.y_bound};
    }

    double ma = a.x_bound.coefficients[1];
    double mb = b.x_bound.coefficients[1];

    double ba = a.x_bound.coefficients[0];
    double bb = b.x_bound.coefficients[0];

    double y = (1/(ma-mb))*(bb - ba);
    double x = ma*y + ba;

    return {x, y};
}

bool valid_triple(Constraint a, Constraint b, Constraint c) {
    if (is_parallel(a, b) && is_parallel(b, c)) {

        auto _a = a.boundary_point();
        auto _b = b.boundary_point();
        auto _c = c.boundary_point();

        if (a.is_lower_bound != b.is_lower_bound) {
            if (!(a.valid_point(_b) && b.valid_point(_a)))
                return false;
        }

        if (a.is_lower_bound != c.is_lower_bound) {
            if (!(a.valid_point(_c) && c.valid_point(_a)))
                return false;
        }

        if (c.is_lower_bound != b.is_lower_bound) {
            if (!(c.valid_point(_b) && b.valid_point(_c)))
                return false;
        }

        return true;
    }

    bool var1 = false;
    bool var2 = false;
    bool var3 = false;

    if (!is_parallel(a, b) && !is_parallel(b, c)) {
        Point ab = intersect(a, b);
        Point bc = intersect(b, c);

        if (approx_zero(ab.dist(bc))) {
            a.shift();
            b.shift();
            c.shift();
        }
    }

    if (!is_parallel(a, b)) {
        auto p = intersect(a, b);
        var1 = c.valid_point(p);
    }

    if (!is_parallel(b, c)) {
        auto p = intersect(b, c);
        var2 = a.valid_point(p);
    }

    if (!is_parallel(a, c)) {
        auto p = intersect(a, c);
        var3 = b.valid_point(p);
    }


    return var1 || var2 || var3;
}

bool valid_constraints(ConstrainedBivariatePolynomial<2>& poly) {
    auto constraints = std::vector<Constraint>();

    if (poly.y_interval.min >= poly.y_interval.max)
        return false;

    constraints.push_back(
        Constraint(true, true, poly.y_interval.min, Polynomial<1>({0, 0}))
    );
    constraints.push_back(
        Constraint(true, false, poly.y_interval.max, Polynomial<1>({0, 0}))
    );

    for (auto f: poly.left_constraints)
        constraints.push_back(
            Constraint(false, true, 0, f)
        );
    
    for (auto f: poly.right_constraints)
        constraints.push_back(
            Constraint(false, false, 0, f)
        );

    for (size_t i = 0; i < constraints.size(); ++i)
        for (size_t j = i+1; j < constraints.size(); ++j) {
            if (constraints[i].zero_between(constraints[j])) {
                return false;
            }
            for (size_t k = j+1; k < constraints.size(); ++k) {
                auto a = constraints[i];
                auto b = constraints[j];        
                auto c = constraints[k];
       
                if ((a!=b && b!=c && c!=a)
                && !valid_triple(a,b,c)) {
                    return false;
                }
            }
        }

    return true;
}

std::vector<Point> sample_domain(ConstrainedBivariatePolynomial<2>& poly) {
    std::vector<Point> output = std::vector<Point>();

    double y_step = (poly.y_interval.max - poly.y_interval.min) / 10;

    for (int i = 0; i < 10; ++i) {
        double y = poly.y_interval.min + i*y_step;
        auto x_interval = poly.interval_at_y(y);
        double x_step = (x_interval.max - x_interval.min) / 10;
        for (int j = 0; j <= 10; ++j) {
            double x = x_interval.min + j*x_step;

            output.push_back({x, y});
        }
    }

    if (poly.right_constraints.size() == 0)
        std::cout << "...\n";

    return output;
}

bool valid_point(ConstrainedBivariatePolynomial<2>& poly, Point p) {
        if (p.y < poly.y_interval.min || p.y > poly.y_interval.max)
            return false;

        for (auto& c: poly.left_constraints) {
            if (p.x < c(p.y))
                return false;
        }

        for (auto& c: poly.right_constraints) {
            if (p.x > c(p.y))
                return false;
        }

        return true;
    }

bool is_positive(ConstrainedBivariatePolynomial<2>& poly) {
    return true;
    auto points = sample_domain(poly);

    for (auto p: points) {
        auto val = poly(p);
        if (valid_point(poly, p) && poly(p) < 0 && !approx_zero(poly(p))) {
            return false;
        }
    }

    return true;
}

double param_space_coord(const Point& s, const Point& t, Point a) {
    Line l = Line::fromTwoPoints(s, t);
    if (!l.includesPoint(a))
        assert(l.includesPoint(a));
    double distance = s.dist(a);
    if (dot(t-s, a-s) < 0)
        distance = -distance;
    return distance;
}

Point param_space(const Point& s1, const Point& t1,
const Point& s2, const Point& t2,
Point a, Point b) {
    return {
        param_space_coord(s1, t1, a),
        param_space_coord(s2, t2, b)
    };
}

bool is_monotone(Line& axis) {
    if (axis.isHorizontal() || axis.isHorizontal())
        return true;
    auto gi = axis.grad_int();
    return gi.x > 0;
}

bool intersect(const Cell& cell, Line v) {
    Point bl = {cell.s.x, cell.s.y};
    Point br = {cell.t.x, cell.s.y};
    Point tl = {cell.s.x, cell.t.y};
    Point tr = {cell.t.x, cell.t.y};
    return v.side(bl) * v.side(tr) == -1 ||
    v.side(tl) * v.side(br) == -1;
}

Line axis_from_point_dir(const Cell& cell, Point l1_p, Point l2_p, Point dir) {
    Point val_p1 = param_space(cell.s1, cell.t1, 
    cell.s2, cell.t2, l1_p, l2_p);
    Point val_p2 = val_p1 + dir;
    return Line::fromTwoPoints(val_p1, val_p2);
}

Line axis_from_points(const Cell& cell, 
Point l1_p1, Point l2_p1, Point l1_p2, Point l2_p2) {
    Point ax1_p1 = param_space(cell.s1, cell.t1, 
    cell.s2, cell.t2, l1_p1, l2_p1);
    Point ax1_p2 = param_space(cell.s1, cell.t1, 
    cell.s2, cell.t2, l1_p2, l2_p2);
    Line axis = Line::fromTwoPoints(ax1_p1, ax1_p2);
    return axis;
}

distance_t l1_h_function(const Cell& cell, Point p) {
    Point vec1 = (cell.t1 - cell.s1) / (cell.t1.dist(cell.s1));
    Point vec2 = (cell.t2 - cell.s2) / (cell.t2.dist(cell.s2));
    Point p1 = cell.s1 + vec1 * p.x;
    Point p2 = cell.s2 + vec2 * p.y;
    return abs(p1.x - p2.x) + abs(p1.y - p2.y);
}


Line get_proper_axis(const Cell& cell, Line v1, Line v2) {
    Point int_point = intersect(v1, v2);
    std::vector<Line> valleys = {v1, v2};
    for (int i = 0; i < 2; ++i) {
       Line v = valleys[i];
       Line other = valleys[(i+1)%2];
       Point upper;
       Point lower;
       if (v.isVertical()) {
            upper = int_point + Point({0, 1});
            lower = int_point + Point({0, -1});
       } else {
            upper = {int_point.x + 1, v.getY(int_point.x + 1)};
            lower = {int_point.x - 1, v.getY(int_point.x - 1)};
       }
        Line upper_line = Line(upper, {1, -1});
        Line lower_line = Line(lower, {1, -1});
        bool upper_valley;
        bool lower_valley;
        Point upper_int;
        if (!isParallel(upper_line, other))
            upper_int = intersect(upper_line, other);
        else
            upper_int = upper + Point({1, -1});
        Point upper_vec = upper_int - upper;
        Point upper_ref = upper - upper_vec;
        upper_valley = l1_h_function(cell, upper_int) >= l1_h_function(cell, upper)
        && l1_h_function(cell, upper_ref) >= l1_h_function(cell, upper);
        
        Point lower_int;
        if (!isParallel(lower_line, other))
            lower_int = intersect(lower_line, other);
        else
            lower_int = lower + Point({1, -1});
        Point lower_vec = lower_int - lower;
        Point lower_ref = lower - lower_vec;
        lower_valley = l1_h_function(cell, lower_int) >= l1_h_function(cell, lower)
        && l1_h_function(cell, lower_ref) >= l1_h_function(cell, lower);
    

        if (upper_valley && lower_valley)
            return v;
    }
    
    throw std::logic_error("neither of the given axes passes the test for a proper axis.");
}


std::vector<Line>
get_axis(const Cell& cell) {
    double sx = cell.s.x; //- cell.mid.x;
    double sy = cell.s.y; //- cell.mid.y;
    double tx = cell.t.x; //- cell.mid.x;
    double ty = cell.t.y; //- cell.mid.y;

    auto result = std::vector<Line>();

    Line l1 = Line::fromTwoPoints(cell.s1, cell.t1);
    Line l2 = Line::fromTwoPoints(cell.s2, cell.t2);

    Point s1 = cell.s1;
    Point s2 = cell.s2;
    Point t1 = cell.t1;
    Point t2 = cell.t2;

    // To determine the valley(s) of the cell, we find a pairs of points
    // in parameter space which define the valleys. For example, the first
    // valley is the line through the pairs (v1_l1_p1, v1_l2_p1) and (v1_l1_p2, v1_l2_p2). 

    // First valley
    Point v1_l1_p1;
    Point v1_l2_p1;
    Point v1_l1_p2;
    Point v1_l2_p2;
    // Second valley
    Point v2_l1_p1;
    Point v2_l2_p1;
    Point v2_l1_p2;
    Point v2_l2_p2;

    bool both_vertical = (l1.isVertical() && l2.isVertical());
    bool equal_slope = approx_equal(l1.grad_int().x, l2.grad_int().x);
    bool parallel = both_vertical || equal_slope;

    if (parallel) {
        if (both_vertical) {
            v1_l1_p1 = {l1.getX(0), 0};
            v1_l2_p1 = {l2.getX(0), 0};
        } else {
            v1_l1_p1 = {0, l1.getY(0)};
            v1_l2_p1 = l2.closest_l1(v1_l1_p1)[0];
        }

        Point dir;
        if (dot(l1.direction, l2.direction) > 0) {
            // In the case of parallel lines, the valley can only have slope 1 ior -1,
            // we don't need a second pair of points. Moreover, we need not consider
            // a valley of slope -1.
            dir = {1, 1};
            Line v = axis_from_point_dir(cell, v1_l1_p1, v1_l2_p1, dir);
            if (intersect(cell, v) && is_monotone(v)) {
                result.push_back(v);
            }
            if (approx_equal(abs(l1.grad_int().x), 1.)) {
                Line v2 = axis_from_point_dir(cell, v1_l1_p1, l2.closest_l1(v1_l1_p1)[1], dir);
                if (intersect(cell, v2) && is_monotone(v2)) {
                    result.push_back(v2);
                }
            }
        } 
    } else {
        Point int_point = intersect(l1, l2);
        Point e1 = {1, 0};
        Point e2 = {0, 1};

        if (l1.isHorizontal()) {
            v1_l1_p1 = int_point-e1;
            v1_l2_p1 = int_point;
            v1_l1_p2 = int_point+e1;
            v1_l2_p2 = int_point;
        } else if (l2.isHorizontal()) {
            v1_l2_p1 = int_point-e1;
            v1_l1_p1 = int_point;
            v1_l2_p2 = int_point+e1;
            v1_l1_p2 = int_point;
        } else {
            if (!l1.isVertical()) {
                v1_l1_p1 = {int_point.x-1, l1.getY(int_point.x-1)};
                v1_l2_p1 = {l2.getX(v1_l1_p1.y), v1_l1_p1.y};
                v1_l1_p2 = {int_point.x+1, l1.getY(int_point.x+1)};
                v1_l2_p2 = {l2.getX(v1_l1_p2.y), v1_l1_p2.y};
            } else if (!l2.isVertical()) {
                v1_l2_p1 = {int_point.x-1, l2.getY(int_point.x-1)};
                v1_l1_p1 = {l1.getX(v1_l2_p1.y), v1_l2_p1.y};
                v1_l2_p2 = {int_point.x+1, l2.getY(int_point.x+1)};
                v1_l1_p2 = {l1.getX(v1_l2_p2.y), v1_l2_p2.y};
            }
        }

        if (l1.isVertical()) {
            v2_l1_p1 = int_point-e2;
            v2_l2_p1 = int_point;
            v2_l1_p2 = int_point+e2;
            v2_l2_p2 = int_point;
        } else if (l2.isVertical()) {
            v2_l2_p1 = int_point-e2;
            v2_l1_p1 = int_point;
            v2_l2_p2 = int_point+e2;
            v2_l1_p2 = int_point;
        } else {
            if (!l1.isHorizontal()) {
                v2_l1_p1 = {l1.getX(int_point.y-1), int_point.y-1};
                v2_l2_p1 = {v2_l1_p1.x, l2.getY(v2_l1_p1.x)};
                v2_l1_p2 = {l1.getX(int_point.y+1), int_point.y+1};
                v2_l2_p2 = {v2_l1_p2.x, l2.getY(v2_l1_p2.x)};
            } else if (!l2.isHorizontal()) {
                v2_l2_p1 = {l2.getX(int_point.y-1), int_point.y-1};
                v2_l1_p1 = {v2_l2_p1.x, l1.getY(v2_l2_p1.x)};
                v2_l2_p2 = {l2.getX(int_point.y+1), int_point.y+1};
                v2_l1_p2 = {v2_l2_p2.x, l1.getY(v2_l2_p2.x)};
            }
        }

        Line v1 = axis_from_points(cell, v1_l1_p1, v1_l2_p1,
            v1_l1_p2, v1_l2_p2);
        Line v2 = axis_from_points(cell, v2_l1_p1, v2_l2_p1,
            v2_l1_p2, v2_l2_p2);

        Line v = get_proper_axis(cell, v1, v2);
        if (intersect(cell, v) && is_monotone(v)) {
            result.push_back(v);
        }
    }
    return result;
}
}