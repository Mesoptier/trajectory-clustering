#ifndef TRAJECTORY_CLUSTERING_2D_L1_L1_H
#define TRAJECTORY_CLUSTERING_2D_L1_L1_H

#include "cdtw.h"
#include "BivariatePolynomial.h"
// #include "1d-l1-l1.h"

//
// 2D + L1 image norm + L1 param norm
//

namespace {

/**
     * Get the angle (measured counter-clockwise) made by the segment with 
     * the horizontal line through the source
    */
    long double angle(Point s, Point t) {
        long double PI = 3.141;

        if (approx_equal(s.x, t.x)) {
            if (s.y > t.y)
                return 3*PI/2;
            else
                return PI/2;
        }

        long double num = s.y - t.y;
        long double denom = s.x - t.x;
        long double arg = fabs(s.y - t.y) / fabs(s.x - t.x);
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

enum VARIABLE_TYPE { X, Y, _ };

struct Constraint {
    bool is_y_constraint;
    bool is_lower_bound;
    double y_bound;
    Polynomial<1> x_bound;

    Constraint(bool y_const, bool lower_bound, double _y_bound, Polynomial<1> _x_bound) :
    is_y_constraint(y_const), is_lower_bound(lower_bound), y_bound(_y_bound), x_bound(_x_bound) {}

    bool valid_point(Point p) {
        if (is_y_constraint) {
            if (is_lower_bound)
                return p.y >= y_bound;
            else
                return p.y <= y_bound;
        } else {
            if (is_lower_bound)
                return p.x >= x_bound.coefficients[1] * p.y + x_bound.coefficients[0];
            else
                return p.x <= x_bound.coefficients[1] * p.y + x_bound.coefficients[0];
        }
    }

    Point normal() {
        if (!is_y_constraint) {

            if (is_lower_bound) {
                return Point(1, -x_bound.coefficients[1]);
            } else {
                return Point(-1, x_bound.coefficients[1]);
            }

        } else {

            if (is_lower_bound) {
                return Point(0, 1);
            } else {
                return Point(0, -1);
            }
        }
    }

    void shift() {

        double eps = 0.01;
        auto n_vec = normal();
        n_vec = Point(eps*n_vec.x, eps*n_vec.y);

        if (!is_y_constraint) {
            x_bound.coefficients = {
                x_bound.coefficients[0] + n_vec.x - x_bound.coefficients[1]*n_vec.y,
                x_bound.coefficients[1]
            };
        } else {
            y_bound += n_vec.y;
        }
    }    

    Point boundary_point() {
        if (is_y_constraint)
            return {0, y_bound};

        if (approx_zero(x_bound.coefficients[1]))
            return {x_bound.coefficients[0], 0};

        return {0, -x_bound.coefficients[0] / x_bound.coefficients[1]};
    }

    bool operator==(const Constraint& rhs) {
        return is_y_constraint == rhs.is_y_constraint
        && is_lower_bound == rhs.is_lower_bound
        && y_bound == rhs.y_bound
        && x_bound == rhs.x_bound;
    }

    bool operator!=(const Constraint& rhs) {
        return !(*this == rhs);
    }
};


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

        // return (b.valid_point(_a) || c.valid_point(_a)) &&
        // (a.valid_point(_b) || c.valid_point(_b)) &&
        // (a.valid_point(_c) || b.valid_point(_c));

        // return (a.is_lower_bound && b.is_lower_bound && c.is_lower_bound)
        // || (!a.is_lower_bound && !b.is_lower_bound && !c.is_lower_bound);
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

        // std::cout << "shift...\n";
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
        for (size_t j = i+1; j < constraints.size(); ++j)
            for (size_t k = j+1; k < constraints.size(); ++k) {
                auto a = constraints[i];
                auto b = constraints[j];        
                auto c = constraints[k];
                // if (k == 5 && j == 4 && i == 2) {
                //     bool valid = valid_triple(a,b,c);
                //     if (!valid)
                //         std::cout << "hello\n";
                // }
                if ((a!=b && b!=c && c!=a)
                && !valid_triple(a,b,c)) {
                    if (approx_equal(poly.f.coefficients[0][1], 0.15683218619384001) && approx_equal(poly.f.coefficients[0][2], 0.70212072617672538))
                        valid_triple(a,b,c);
                    // std::cout << k << " " << j << " " << i << std::endl;
                    return false;
                }
            }

    return true;
}

/**
 * 
 * Samples points from the domian
 * of a constrained bivariate quadratic
 */
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
            // if (approx_equal(poly.f.coefficients[0][0], -0.001356284712799239) && approx_equal(poly.f.coefficients[0][1], 0.021595794862823228))
            //     std::cout << "hi\n";
            return false;
        }
    }

    return true;
}

/**
 * Explicity solves an integral and returns a vector
 * of constrained bivariate quadratic functions
 * 
 * 
*/
std::vector<ConstrainedBivariatePolynomial<2>>
solve_integral(BivariatePolynomial<1> low_lim, 
BivariatePolynomial<1> hi_lim,
BivariatePolynomial<1> integrand, double t_coeff,
Interval y_range, std::vector<Polynomial<1>> left_constraints, std::vector<Polynomial<1>> right_constraints) {
    auto output = std::vector<ConstrainedBivariatePolynomial<2>>();

        // if (approx_equal(y_range.max, 0.18449462693244933)) {
        //     std::cout << "hi\n";
        // }
        
        /**
         * Impose low_lim <= hi_lim constraints
        */
        double lo_a = low_lim.coefficients[1][0];
        double lo_b = low_lim.coefficients[0][1];
        double lo_c = low_lim.coefficients[0][0];

        double hi_a = hi_lim.coefficients[1][0];
        double hi_b = hi_lim.coefficients[0][1];
        double hi_c = hi_lim.coefficients[0][0];


        if (hi_a > lo_a) {
            left_constraints.push_back(Polynomial<1>({
                (lo_c - hi_c) / (hi_a - lo_a), (lo_b - hi_b) / (hi_a - lo_a) 
            }));
        } else if (hi_a < lo_a) {
            right_constraints.push_back(Polynomial<1>({
                (lo_c - hi_c) / (hi_a - lo_a), (lo_b - hi_b) / (hi_a - lo_a) 
            }));
        } else {
            if (hi_b > lo_b) {
                y_range.min = std::max(y_range.min, (lo_c - hi_c) / (hi_b - lo_b));
            } else if (hi_b < lo_b) {
                y_range.max = std::min(y_range.max, (lo_c - hi_c) / (hi_b - lo_b));
            } else {
                if (hi_c < lo_c) {
                    y_range.max = -1;
                }
            }
        }

        auto p1x = low_lim.coefficients[1][0];
        auto p1y = low_lim.coefficients[0][1];
        auto p1c = low_lim.coefficients[0][0];

        auto p2x = hi_lim.coefficients[1][0];
        auto p2y = hi_lim.coefficients[0][1];
        auto p2c = hi_lim.coefficients[0][0];

        auto p3x = integrand.coefficients[1][0];
        auto p3y = integrand.coefficients[0][1];
        auto p3c = integrand.coefficients[0][0];

        auto a = t_coeff;

        //case 1
        auto case_1_poly = BivariatePolynomial<2>({{{{a*pow(p1c, 2)/2 - a*pow(p2c, 2)/2 - p1c*p3c + p2c*p3c,a*p1c*p1y - a*p2c*p2y - p1c*p3y - p1y*p3c + p2c*p3y + p2y*p3c,a*pow(p1y, 2)/2 - a*pow(p2y, 2)/2 - p1y*p3y + p2y*p3y}},
{{a*p1c*p1x - a*p2c*p2x - p1c*p3x - p1x*p3c + p2c*p3x + p2x*p3c,a*p1x*p1y - a*p2x*p2y - p1x*p3y - p1y*p3x + p2x*p3y + p2y*p3x,0}},
{{a*pow(p1x, 2)/2 - a*pow(p2x, 2)/2 - p1x*p3x + p2x*p3x,0,0}}}});
        
        auto case_1_right_constraints = std::vector<Polynomial<1>>();
        auto case_1_left_constraints = std::vector<Polynomial<1>>();

        double y_upper_bound = y_range.max;
        double y_lower_bound = y_range.min;

        if (p3x > a*p2x) {
            case_1_left_constraints.push_back(Polynomial<1>({
                (a*p2c - p3c) / (p3x-a*p2x), (a*p2y-p3y) / (p3x-a*p2x) 
            }));
        } else if (p3x < a*p2x) {
            case_1_right_constraints.push_back(Polynomial<1>({
                (a*p2c-p3c) / (p3x-a*p2x), (a*p2y-p3y) / (p3x-a*p2x) 
            }));
        } else if (approx_equal(p3x, a*p2x)) {
            if (a*p2y > p3y)
                y_upper_bound = std::min(y_upper_bound, (p3c-a*p2c) / (a*p2y-p3y));
            else if (a*p2y < p3y)
                y_lower_bound = std::max(y_lower_bound, (p3c-a*p2c) / (a*p2y-p3y));
            else if (approx_equal(a*p2y, p3y)) {
                if (!(p3c-a*p2c >= 0)) {
                    
                    y_upper_bound = -1;
                }
            }
        }

        if (p3x > a*p1x) {
            case_1_left_constraints.push_back(Polynomial<1>({
                (a*p1c-p3c) / (p3x-a*p1x), (a*p1y-p3y) / (p3x-a*p1x) 
            }));
        } else if (p3x < a*p1x) {
            case_1_right_constraints.push_back(Polynomial<1>({
                (a*p1c-p3c) / (p3x-a*p1x), (a*p1y-p3y) / (p3x-a*p1x) 
            }));
        } else if (approx_equal(p3x, a*p1x)) {
            if (a*p1y > p3y)
                y_upper_bound = std::min(y_upper_bound, (p3c-a*p1c) / (a*p1y-p3y));
            else if (a*p1y < p3y)
                y_lower_bound = std::max(y_lower_bound, (p3c-a*p1c) / (a*p1y-p3y));
            else if (approx_equal(a*p1y, p3y)) {
                if (!(p3c-a*p1c >= 0)) {
                    y_upper_bound = -1;
                }
            }
        }

        for (auto constraint: left_constraints) 
            case_1_left_constraints.push_back(constraint);

        for (auto constraint: right_constraints) 
            case_1_right_constraints.push_back(constraint);
        
        output.push_back(ConstrainedBivariatePolynomial<2>{
            case_1_poly,
            {y_lower_bound, y_upper_bound},
            case_1_left_constraints,
            case_1_right_constraints
        });

        //case 2
        auto case_2_poly = BivariatePolynomial<2>({{{{-a*pow(p1c, 2)/2 + a*pow(p2c, 2)/2 + p1c*p3c - p2c*p3c,-a*p1c*p1y + a*p2c*p2y + p1c*p3y + p1y*p3c - p2c*p3y - p2y*p3c,-a*pow(p1y, 2)/2 + a*pow(p2y, 2)/2 + p1y*p3y - p2y*p3y}},
{{-a*p1c*p1x + a*p2c*p2x + p1c*p3x + p1x*p3c - p2c*p3x - p2x*p3c,-a*p1x*p1y + a*p2x*p2y + p1x*p3y + p1y*p3x - p2x*p3y - p2y*p3x,0}},
{{-a*pow(p1x, 2)/2 + a*pow(p2x, 2)/2 + p1x*p3x - p2x*p3x,0,0}}}});
        
        auto case_2_right_constraints = std::vector<Polynomial<1>>();
        auto case_2_left_constraints = std::vector<Polynomial<1>>();

        y_upper_bound = y_range.max;
        y_lower_bound = y_range.min;

        if (p3x > a*p2x) {
            case_2_right_constraints.push_back(Polynomial<1>({
                (a*p2c-p3c) / (p3x-a*p2x), (a*p2y-p3y) / (p3x-a*p2x) 
            }));
        } else if (p3x < a*p2x) {
            case_2_left_constraints.push_back(Polynomial<1>({
                (a*p2c-p3c) / (p3x-a*p2x), (a*p2y-p3y) / (p3x-a*p2x) 
            }));
        } else if (approx_equal(p3x, a*p2x)) {
            if (a*p2y > p3y)
                y_lower_bound = std::max(y_lower_bound, (a*p2c-p3c) / (p3y-a*p2y));
            else if (a*p2y < p3y)
                y_upper_bound = std::min(y_upper_bound, (a*p2c-p3c) / (p3y-a*p2y));
            else if (approx_equal(a*p2y, p3y)) {
                if (!(p3c-a*p2c <= 0)) {
                    y_upper_bound = -1;
                }
            }
        }

        if (p3x > a*p1x) {
            case_2_right_constraints.push_back(Polynomial<1>({
                (a*p1c-p3c) / (p3x-a*p1x), (a*p1y-p3y) / (p3x-a*p1x) 
            }));
        } else if (p3x < a*p1x) {
            case_2_left_constraints.push_back(Polynomial<1>({
                (a*p1c-p3c) / (p3x-a*p1x), (a*p1y-p3y) / (p3x-a*p1x) 
            }));
        } else if (approx_equal(p3x, a*p1x)) {
            if (a*p1y > p3y)
                y_lower_bound = std::max(y_lower_bound, (p3c-a*p1c) / (a*p1y-p3y));
            else if (a*p1y < p3y)
                y_upper_bound = std::min(y_upper_bound, (p3c-a*p1c) / (a*p1y-p3y));
            else if (approx_equal(a*p1y, p3y)) {
                if (!(p3c-a*p1c <= 0)) {
                    y_upper_bound = -1;
                }
            }
        }

        for (auto constraint: left_constraints) 
            case_2_left_constraints.push_back(constraint);

        for (auto constraint: right_constraints) 
            case_2_right_constraints.push_back(constraint);
        
        output.push_back(ConstrainedBivariatePolynomial<2>{
            case_2_poly,
            {y_lower_bound, y_upper_bound},
            case_2_left_constraints,
            case_2_right_constraints
        });
        
        //case 3a
        auto case_3a_poly = BivariatePolynomial<2>(
{{{{(pow(a, 2)*pow(p1c, 2) + pow(a, 2)*pow(p2c, 2) - 2*a*p1c*p3c - 2*a*p2c*p3c + 2*pow(p3c, 2))/(2*a),(pow(a, 2)*p1c*p1y + pow(a, 2)*p2c*p2y - a*p1c*p3y - a*p1y*p3c - a*p2c*p3y - a*p2y*p3c + 2*p3c*p3y)/a,(pow(a, 2)*pow(p1y, 2) + pow(a, 2)*pow(p2y, 2) - 2*a*p1y*p3y - 2*a*p2y*p3y + 2*pow(p3y, 2))/(2*a)}},
{{(pow(a, 2)*p1c*p1x + pow(a, 2)*p2c*p2x - a*p1c*p3x - a*p1x*p3c - a*p2c*p3x - a*p2x*p3c + 2*p3c*p3x)/a,(pow(a, 2)*p1x*p1y + pow(a, 2)*p2x*p2y - a*p1x*p3y - a*p1y*p3x - a*p2x*p3y - a*p2y*p3x + 2*p3x*p3y)/a,0}},
{{(pow(a, 2)*pow(p1x, 2) + pow(a, 2)*pow(p2x, 2) - 2*a*p1x*p3x - 2*a*p2x*p3x + 2*pow(p3x, 2))/(2*a),0,0}}}});
        
        auto case_3a_right_constraints = std::vector<Polynomial<1>>();
        auto case_3a_left_constraints = std::vector<Polynomial<1>>();

        y_upper_bound = y_range.max;
        y_lower_bound = y_range.min;

        if (p3x > a*p2x) {
            case_3a_right_constraints.push_back(Polynomial<1>({
                (a*p2c-p3c) / (p3x-a*p2x), (a*p2y-p3y) / (p3x-a*p2x) 
            }));
        } else if (p3x < a*p2x) {
            case_3a_left_constraints.push_back(Polynomial<1>({
                (a*p2c-p3c) / (p3x-a*p2x), (a*p2y-p3y) / (p3x-a*p2x) 
            }));
        } else if (approx_equal(p3x, p2x)) {
            if (a*p2y > p3y)
                y_lower_bound = std::max(y_lower_bound, (p3c-a*p2c) / (a*p2y-p3y));
            else if (a*p2y < p3y)
                y_upper_bound = std::min(y_upper_bound, (p3c-a*p2c) / (a*p2y-p3y));
            else if (approx_equal(a*p2y, p3y)) {
                if (!(p3c-a*p2c <= 0)) {
                    y_upper_bound = -1;
                }
            }
        }

        if (p3x > a*p1x) {
            case_3a_left_constraints.push_back(Polynomial<1>({
                (a*p1c-p3c) / (p3x-a*p1x), (a*p1y-p3y) / (p3x-a*p1x) 
            }));
        } else if (p3x < a*p1x) {
            case_3a_right_constraints.push_back(Polynomial<1>({
                (a*p1c-p3c) / (p3x-a*p1x), (a*p1y-p3y) / (p3x-a*p1x) 
            }));
        } else if (approx_equal(p3x, a*p1x)) {
            if (a*p1y > p3y)
                y_upper_bound = std::min(y_upper_bound, (p3c-a*p1c) / (a*p1y-p3y));
            else if (a*p1y < p3y)
                y_lower_bound = std::max(y_lower_bound, (p3c-a*p1c) / (a*p1y-p3y));
            else if (approx_equal(a*p1y, p3y)) {
                if (!(p3c-a*p1c >= 0)) {
                    y_upper_bound = -1;
                }
            }
        }

        for (auto constraint: left_constraints) 
            case_3a_left_constraints.push_back(constraint);

         for (auto constraint: right_constraints) 
            case_3a_right_constraints.push_back(constraint);

        output.push_back(ConstrainedBivariatePolynomial<2>{
            case_3a_poly,
            {y_lower_bound, y_upper_bound},
            case_3a_left_constraints,
            case_3a_right_constraints
        });

        //case 3b
        auto case_3b_poly = BivariatePolynomial<2>({{{{(-pow(a, 2)*pow(p1c, 2) - pow(a, 2)*pow(p2c, 2) + 2*a*p1c*p3c + 2*a*p2c*p3c - 2*pow(p3c, 2))/(2*a),(-pow(a, 2)*p1c*p1y - pow(a, 2)*p2c*p2y + a*p1c*p3y + a*p1y*p3c + a*p2c*p3y + a*p2y*p3c - 2*p3c*p3y)/a,(-pow(a, 2)*pow(p1y, 2) - pow(a, 2)*pow(p2y, 2) + 2*a*p1y*p3y + 2*a*p2y*p3y - 2*pow(p3y, 2))/(2*a)}},
{{(-pow(a, 2)*p1c*p1x - pow(a, 2)*p2c*p2x + a*p1c*p3x + a*p1x*p3c + a*p2c*p3x + a*p2x*p3c - 2*p3c*p3x)/a,(-pow(a, 2)*p1x*p1y - pow(a, 2)*p2x*p2y + a*p1x*p3y + a*p1y*p3x + a*p2x*p3y + a*p2y*p3x - 2*p3x*p3y)/a,0}},
{{(-pow(a, 2)*pow(p1x, 2) - pow(a, 2)*pow(p2x, 2) + 2*a*p1x*p3x + 2*a*p2x*p3x - 2*pow(p3x, 2))/(2*a),0,0}}}});

        auto case_3b_right_constraints = std::vector<Polynomial<1>>();
        auto case_3b_left_constraints = std::vector<Polynomial<1>>();

        y_upper_bound = y_range.max;
        y_lower_bound = y_range.min;

        if (p3x > a*p2x) {
            case_3b_left_constraints.push_back(Polynomial<1>({
                (a*p2c-p3c) / (p3x-a*p2x), (a*p2y-p3y) / (p3x-a*p2x) 
            }));
        } else if (p3x < p2x) {
            case_3b_right_constraints.push_back(Polynomial<1>({
                (a*p2c-p3c) / (p3x-a*p2x), (a*p2y-p3y) / (p3x-a*p2x) 
            }));
        } else if (approx_equal(p3x, a*p2x)) {
            if (a*p2y > p3y)
                y_upper_bound = std::min(y_upper_bound, (p3c-a*p2c) / (a*p2y-p3y));
            else if (a*p2y < p3y)
                y_lower_bound = std::max(y_lower_bound, (p3c-a*p2c) / (a*p2y-p3y));
            else if (approx_equal(a*p2y, p3y)) {
                if (!(p3c-a*p2c >= 0)) {
                    y_upper_bound = -1;
                }
            }    
        }

        if (p3x > a*p1x) {
            case_3b_right_constraints.push_back(Polynomial<1>({
                (a*p1c-p3c) / (p3x-a*p1x), (a*p1y-p3y) / (p3x-a*p1x) 
            }));
        } else if (p3x < a*p1x) {
            case_3b_left_constraints.push_back(Polynomial<1>({
                (a*p1c-p3c) / (p3x-a*p1x), (a*p1y-p3y) / (p3x-a*p1x) 
            }));
        } else if (approx_equal(p3x, a*p1x)) {
            if (a*p1y > p3y)
                y_lower_bound = std::max(y_lower_bound, (p3c-a*p1c) / (a*p1y-p3y));
            else if (a*p1y < p3y)
                y_upper_bound = std::min(y_upper_bound, (p3c-a*p1c) / (a*p1y-p3y));
            else if (approx_equal(a*p1y, p3y)) {
                if (!(p3c-a*p1c <= 0)) {
                    y_upper_bound = -1;
                }
            }
        }

        for (auto constraint: left_constraints) 
            case_3b_left_constraints.push_back(constraint);

         for (auto constraint: right_constraints) 
            case_3b_right_constraints.push_back(constraint);

        output.push_back(ConstrainedBivariatePolynomial<2>{
            case_3b_poly,
            {y_lower_bound, y_upper_bound},
            case_3b_left_constraints,
            case_3b_right_constraints
        });

    auto final_output = std::vector<ConstrainedBivariatePolynomial<2>>();

        // if (approx_equal(y_range.max, 0.18449462693244933)) {
        //     std::cout << "hi\n";
        // }
        

    for (auto& poly: output)
        if (valid_constraints(poly)) {
            if (!is_positive(poly)) {
                std::cout << "negative\n";
                auto res = is_positive(poly);
                valid_constraints(poly);
            } 
            poly.clean_constraints();
            final_output.push_back(poly);
        }

    
    return final_output;
}


//TODO: Implement this properly
bool has_positive_slope(Line line) {
    return false;
}

distance_t get_angle(Line line) {
    return 0;
}

bool contains_midpoint(const Cell& cell) {
    // auto integral = solve_integral(
    // X, 1, 1, 1, 1, Y,
    // 1, Interval{0, 1}, Interval{0, 1}
    // );
    auto& len1 = cell.len1;
    auto& len2 = cell.len2;
    auto& s1 = cell.s1;
    auto& s2 = cell.s2;
    auto& t1 = cell.t1;
    auto& t2 = cell.t2;

    const auto l1 = Line::fromTwoPoints(s1, t1);
    const auto l2 = Line::fromTwoPoints(s2, t2);

    if (isParallel(l1, l2))
        return false;

    auto& mid = cell.mid;

    return (mid.x*(-len1 + mid.x) < 0) 
    && (mid.y*(-len2 + mid.y) < 0)
    && !approx_equal(s1, s2)
    && !approx_equal(t1, t2);    
}

std::vector<Cell> 
get_subdivision(const Cell& cell) {
    
    auto result = std::vector<Cell>();

    if (!contains_midpoint(cell))
        return result;

    auto& mid = cell.mid;

    Cell c1 = Cell(cell.s1, cell.s2, mid, mid);
    Cell c2 = Cell(cell.s1, mid, mid, cell.t2);
    Cell c3 = Cell(mid, cell.s2, cell.t1, mid);
    Cell c4 = Cell(mid, mid, cell.t1, cell.t2);

    result.push_back(c1);
    result.push_back(c2);
    result.push_back(c3);
    result.push_back(c4);

    return result;
}

double param_space_coord(const Point& s, const Point& t, Point a) {
    Line l = Line::fromTwoPoints(s, t);

    if (!l.includesPoint(a))
        assert(l.includesPoint(a));

    double distance = s.dist(a);

    if ((s.x - t.x) * (s.x - a.x) >= 0)
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

/**
 * Determines if axis is monote
 * @param axis
 * @return true if axis is monotone, false otherwise
 * */
bool is_monotone(Line& axis) {
    auto gi = axis.grad_int();
    return gi.x > 0;
}

struct CellIntersections {
    Point left;
    Point right;
    Point bottom;
    Point top;

    std::vector<bool> int_exists;
};

/**
 * computes the intersections
 * of the given line with the four lines
 * through the cell boundaries
 */
CellIntersections get_intersections(const Cell& cell, Line& line) {
    double sx = cell.s.x; //- cell.mid.x;
    double sy = cell.s.y; //- cell.mid.y;
    double tx = cell.t.x; //- cell.mid.x;
    double ty = cell.t.y; //- cell.mid.y;

    auto left_line = Line::fromTwoPoints({sx, sy}, {sx, ty});
    auto right_line = Line::fromTwoPoints({tx, sy}, {tx, ty});
    auto bottom_line = Line::fromTwoPoints({sx, sy}, {tx, sy});
    auto top_line = Line::fromTwoPoints({sx, ty}, {tx, ty});

    CellIntersections result = CellIntersections();
    result.int_exists = std::vector<bool>(4);


    if (!isParallel(line, left_line)) {
        result.int_exists[0] = true;
        result.left = intersect(line, left_line);
    } else {
        result.int_exists[0] = false;
    }

    if (!isParallel(line, right_line)) {
        result.int_exists[1] = true;
        result.right = intersect(line, right_line);
    } else {
        result.int_exists[1] = false;
    }

    if (!isParallel(line, bottom_line)) {
        result.int_exists[2] = true;
        result.bottom = intersect(line, bottom_line);
    } else {
        result.int_exists[2] = false;
    }

    if (!isParallel(line, top_line)) {
        result.int_exists[3] = true;
        result.top = intersect(line, top_line);
    } else {
        result.int_exists[3] = false;
    }

    return result;
}

bool intersects_cell(const Cell& cell, CellIntersections& intersections) {
    double sx = cell.s.x; //- cell.mid.x;
    double sy = cell.s.y; //- cell.mid.y;
    double tx = cell.t.x; //- cell.mid.x;
    double ty = cell.t.y; //- cell.mid.y;

    int int_count = 0;

    bool left_int = intersections.int_exists[0]
    && (intersections.left.y >= sy && intersections.left.y <= ty);

    if (left_int)
        int_count++;

    bool right_int = intersections.int_exists[1]
    && (intersections.right.y >= sy && intersections.right.y <= ty);

    if (right_int)
        int_count++;
    
    bool bottom_int = intersections.int_exists[2]
    && (intersections.bottom.x >= sx && intersections.bottom.x <= tx);


    if (bottom_int)
        int_count++;
    
    bool top_int = intersections.int_exists[3]
    && (intersections.top.x >= sx && intersections.top.y <= tx);

    if (top_int)
        int_count++;

    return int_count >= 2;
}

/**
* Get axes for given cell
* @param cell
* @return vector of lines corresponding to the axes which
* intersect the cell
*/
std::vector<Line>
get_axes(const Cell& cell) {
    // std::cout << "getting axes...\n";
    double sx = cell.s.x; //- cell.mid.x;
    double sy = cell.s.y; //- cell.mid.y;
    double tx = cell.t.x; //- cell.mid.x;
    double ty = cell.t.y; //- cell.mid.y;

    auto result = std::vector<Line>();

    Line l1 = Line::fromTwoPoints(cell.s1, cell.t1);
    Line l2 = Line::fromTwoPoints(cell.s2, cell.t2);

    Point gi_1 = l1.grad_int();
    Point gi_2 = l2.grad_int();

    double m1 = gi_1.x;
    double m2 = gi_2.x;
    double b1 = gi_1.y;
    double b2 = gi_2.y;

    // std::cout << m1 << std::endl;
    // std::cout << m2 << std::endl;
    // std::cout << b1 << std::endl;
    // std::cout << b2 << std::endl;

    double ax1_y1 = (b1/m1 - b2/m2 + 1) / (1/m1 - 1/m2);
    double ax1_y2 = (b1/m1 - b2/m2 - 1) / (1/m1 - 1/m2);

    double ax1_l1_x1 = l1.getX(ax1_y1);
    double ax1_l2_x1 = l2.getX(ax1_y1);
    double ax1_l1_x2 = l1.getX(ax1_y2);
    double ax1_l2_x2 = l2.getX(ax1_y2);

    Point ax1_p1 = param_space(cell.s1, cell.t1, cell.s2, cell.t2,
    {ax1_l1_x1, ax1_y1}, {ax1_l2_x1, ax1_y1});
    Point ax1_p2 = param_space(cell.s1, cell.t1, cell.s2, cell.t2,
    {ax1_l1_x2, ax1_y2}, {ax1_l2_x2, ax1_y2});

    Line axis_1 = Line::fromTwoPoints(ax1_p1, ax1_p2);
    // if (!axis_1.includesPoint(cell.mid))
    //     assert(false);

    auto axis_1_intersections = get_intersections(cell, axis_1);

    if (intersects_cell(cell, axis_1_intersections) && is_monotone(axis_1)) {
        auto translated_axis_1 = Line::fromTwoPoints(ax1_p1 - cell.mid , ax1_p2 - cell.mid);
        result.push_back(translated_axis_1);
    }

    double ax2_x1 = (b2 - b1 + 1) / (m1 - m2);
    double ax2_x2 = (b2 - b1 - 1) / (m1 - m2);

    double ax2_l1_y1 = l1.getY(ax2_x1);
    double ax2_l2_y1 = l2.getY(ax2_x1);
    double ax2_l1_y2 = l1.getY(ax2_x2);
    double ax2_l2_y2 = l2.getY(ax2_x2);

    Point ax2_p1 = param_space(cell.s1, cell.t1, cell.s2, cell.t2,
    {ax2_x1, ax2_l1_y1}, {ax2_x1, ax2_l2_y1});
    Point ax2_p2 = param_space(cell.s1, cell.t1, cell.s2, cell.t2,
    {ax2_x2, ax2_l1_y2}, {ax2_x2, ax2_l2_y2});

    Line axis_2 = Line::fromTwoPoints(ax2_p1, ax2_p2);
    // if (!axis_2.includesPoint(cell.mid))
    //     assert(false);

    auto axis_2_intersections = get_intersections(cell, axis_2);

    if (intersects_cell(cell, axis_2_intersections) && is_monotone(axis_2)) {
        auto translated_axis_2 = Line::fromTwoPoints(ax2_p1 - cell.mid , ax2_p2 - cell.mid);
        result.push_back(translated_axis_2);
    }

    // for (auto l: result) {
    //     if (!l.includesPoint({0, 0}))
    //         assert(l.includesPoint({0, 0}));
    // }

    // std::cout << "intersecting axes: " << result.size() << std::endl;

    return result;
}

std::vector<ConstrainedBivariatePolynomial<2>>
horizontal_int(BivariatePolynomial<1> low_lim, BivariatePolynomial<1> hi_lim,
const Cell& cell, double y_coef, bool y_const) {
    std::vector<ConstrainedBivariatePolynomial<2>> xterms;
    std::vector<ConstrainedBivariatePolynomial<2>> yterms;

    std::vector<ConstrainedBivariatePolynomial<2>> results;

    Point s1 = cell.s1;
    Point s2 = cell.s2;
    Point t1 = cell.t1;
    Point t2 = cell.t2;

    double sx = cell.s.x; //- cell.mid.x;
    double sy = cell.s.y; //- cell.mid.y;
    double tx = cell.t.x; //- cell.mid.x;
    double ty = cell.t.y; //- cell.mid.y;

    double alpha = angle(s1, t1);
    double beta = angle(s2, t2);

    if (y_const) {
        auto terms_x = solve_integral(low_lim, hi_lim,
            BivariatePolynomial<1>({{{{s1.x-s2.x - y_coef*cos(beta), 0}},{{0, 0}}}}),
            cos(alpha), {sy, ty},
            {Polynomial<1>({sx, 0})}, {Polynomial<1>({tx, 0})} 
        );

        for (auto term: terms_x)  
            xterms.push_back(term);

        auto terms_y = solve_integral(low_lim, hi_lim,
            BivariatePolynomial<1>({{{{s1.y-s2.y - y_coef*sin(beta), 0}},{{0, 0}}}}),
            sin(alpha), {sy, ty},
            {Polynomial<1>({sx, 0})}, {Polynomial<1>({tx, 0})} 
        );

        for (auto term: terms_y)  
            yterms.push_back(term);

            
    } else {
        auto terms_x = solve_integral(low_lim, hi_lim,
            BivariatePolynomial<1>({{{{s1.x-s2.x, -cos(beta)}},{{0, 0}}}}),
            cos(alpha), {sy, ty},
            {Polynomial<1>({sx, 0})}, {Polynomial<1>({tx, 0})} 
        );

        for (auto term: terms_x)  
            xterms.push_back(term);

        auto terms_y = solve_integral(low_lim, hi_lim,
            BivariatePolynomial<1>({{{{s1.y-s2.y, -sin(beta)}},{{0, 0}}}}),
            sin(alpha), {sy, ty},
            {Polynomial<1>({sx, 0})}, {Polynomial<1>({tx, 0})} 
        );

        for (auto term: terms_y)  
            yterms.push_back(term);
        
    }

    for (auto xterm: xterms)
        for (auto yterm: yterms)
            results.push_back(xterm + yterm);

    auto valid_results = std::vector<ConstrainedBivariatePolynomial<2>>();

    for (auto f: results)
        if (valid_constraints(f))
            valid_results.push_back(f);

    return valid_results;
}

std::vector<ConstrainedBivariatePolynomial<2>>
vertical_int(BivariatePolynomial<1> low_lim, BivariatePolynomial<1> hi_lim,
const Cell& cell, double x_coef, bool x_const
) {

    std::vector<ConstrainedBivariatePolynomial<2>> xterms;
    std::vector<ConstrainedBivariatePolynomial<2>> yterms;

    std::vector<ConstrainedBivariatePolynomial<2>> results;

    Point s1 = cell.s1;
    Point s2 = cell.s2;
    Point t1 = cell.t1;
    Point t2 = cell.t2;

    double sx = cell.s.x; //- cell.mid.x;
    double sy = cell.s.y; //- cell.mid.y;
    double tx = cell.t.x; //- cell.mid.x;
    double ty = cell.t.y; //- cell.mid.y;

    double alpha = angle(s1, t1);
    double beta = angle(s2, t2);

    
    if (x_const) {
        auto terms_x = solve_integral(low_lim, hi_lim,
            BivariatePolynomial<1>({{{{x_coef*cos(alpha) + s1.x-s2.x, 0}},{{0, 0}}}}),
            cos(beta), {sy, ty},
            {Polynomial<1>({sx, 0})}, {Polynomial<1>({tx, 0})}
        );

        for (auto term: terms_x)
            xterms.push_back(term);


        auto terms_y = solve_integral(low_lim, hi_lim,
            BivariatePolynomial<1>({{{{x_coef * sin(alpha) + s1.y-s2.y, 0}},{{0, 0}}}}),
            sin(beta), {sy, ty},
            {Polynomial<1>({sx, 0})}, {Polynomial<1>({tx, 0})}
        );

        for (auto term: terms_y)
            yterms.push_back(term);
    
    } else {
    
        auto terms_x = solve_integral(low_lim, hi_lim,
            BivariatePolynomial<1>({{{{s1.x-s2.x, 0}},{{cos(alpha), 0}}}}),
            cos(beta), {sy, ty},
            {Polynomial<1>({sx, 0})}, {Polynomial<1>({tx, 0})}
        );

        for (auto term: terms_x)
            xterms.push_back(term);
            
        auto terms_y = solve_integral(low_lim, hi_lim,
            BivariatePolynomial<1>({{{{s1.y-s2.y, 0}},{{sin(alpha), 0}}}}),
            sin(beta), {sy, ty},
            {Polynomial<1>({sx, 0})}, {Polynomial<1>({tx, 0})}
        );

        for (auto term: terms_y)
            yterms.push_back(term);
    }

    for (auto xterm: xterms)
        for (auto yterm: yterms)
            results.push_back(xterm + yterm);

    auto valid_results = std::vector<ConstrainedBivariatePolynomial<2>>();

    for (auto f: results)
        if (valid_constraints(f))
            valid_results.push_back(f);

    return valid_results;
}


std::vector<ConstrainedBivariatePolynomial<2>>
axis_int(BivariatePolynomial<1> low_lim, BivariatePolynomial<1> hi_lim, const Cell& cell, Line axis,
std::vector<Polynomial<1>> left_constraints, std::vector<Polynomial<1>> right_constraints, 
Interval y_range) {
    Point s1 = cell.s1;
    Point s2 = cell.s2;
    Point t1 = cell.t1;
    Point t2 = cell.t2;

    double sx = cell.s.x; //- cell.mid.x;
    double sy = cell.s.y; //- cell.mid.y;
    double tx = cell.t.x; //- cell.mid.x;
    double ty = cell.t.y; //- cell.mid.y;

    double alpha = angle(s1, t1);
    double beta = angle(s2, t2);

    auto gi = axis.grad_int();

    double m = gi.x;
    double b = gi.y;

    /*
    * int_{l0}^{l1} | h(t, mt+b) |dt = int | t*cos(alpha) - (mt+b)*cos(beta) | + |t*sin(alpha) - (mt+b)*sin(alpha)|dt
    */

    auto xterms = solve_integral(
        low_lim, hi_lim,
        BivariatePolynomial<1>({{{{s1.x-s2.x - b*cos(beta), 0}},{{0, 0}}}}),
        m*cos(beta) - cos(alpha), y_range,
        left_constraints, right_constraints
    );

    auto yterms = solve_integral(
        low_lim, hi_lim,
        BivariatePolynomial<1>({{{{s1.y-s2.y - b*sin(beta), 0}},{{0, 0}}}}),
        m*sin(beta) - sin(alpha), y_range,
        left_constraints, right_constraints
    );

    std::vector<ConstrainedBivariatePolynomial<2>> results = std::vector<ConstrainedBivariatePolynomial<2>>();

    for (auto xterm: xterms)
        for (auto yterm: yterms)
            results.push_back(xterm + yterm);

    auto valid_results = std::vector<ConstrainedBivariatePolynomial<2>>();

    for (auto f: results)
        if (valid_constraints(f)) {
            if (!is_positive(f))
                std::cout << "why\n";
            valid_results.push_back(f.multiply(fabs(1+m)));
        }

    return valid_results;
}

std::vector<ConstrainedBivariatePolynomial<2>>
combine_steps(size_t step_count, std::vector<std::vector<ConstrainedBivariatePolynomial<2>>> steps) {
    
    auto results = std::vector<ConstrainedBivariatePolynomial<2>>();
    
    if (step_count == 2) 
        for (size_t i = 0; i < steps[0].size(); ++i)
            for (size_t j = 0; j < steps[1].size(); ++j)
                results.push_back(steps[0][i] + steps[1][j]);
    else if (step_count == 3)
        for (size_t i = 0; i < steps[0].size(); ++i)
            for (size_t j = 0; j < steps[1].size(); ++j)
                for (size_t k = 0; k < steps[2].size(); ++k) {
                    auto lhs = steps[0][i] + steps[1][j];
                    results.push_back(lhs + steps[2][k]);
                }


    auto valid_results = std::vector<ConstrainedBivariatePolynomial<2>>();

    for (auto f: results)
        if (valid_constraints(f)) {
            f.clean_constraints();
            valid_results.push_back(f);
        }

    return valid_results;
}

std::vector<ConstrainedBivariatePolynomial<2>>
combine_terms(size_t step_count, 
std::vector<std::vector<ConstrainedBivariatePolynomial<2>>> xterms,
std::vector<std::vector<ConstrainedBivariatePolynomial<2>>> yterms) {
    assert(step_count == xterms.size());
    assert(step_count == yterms.size());

    auto results = std::vector<ConstrainedBivariatePolynomial<2>>();

    auto combined_steps = std::vector<std::vector<ConstrainedBivariatePolynomial<2>>>(step_count);

    
    for (size_t i = 0; i < step_count; ++i)
        for (size_t j = 0; j < xterms[i].size(); ++j)
            for (size_t k = 0; k < yterms[i].size(); ++k)
                combined_steps[i].push_back(xterms[i][j]+yterms[i][k]);

    if (step_count == 2) 
        for (size_t i = 0; i < combined_steps[0].size(); ++i)
            for (size_t j = 0; j < combined_steps[1].size(); ++j)
                results.push_back(combined_steps[0][i] + combined_steps[1][j]);
    else if (step_count == 3)
        for (size_t i = 0; i < combined_steps[0].size(); ++i)
            for (size_t j = 0; j < combined_steps[1].size(); ++j)
                for (size_t k = 0; k < combined_steps[2].size(); ++k) {
                    auto lhs = combined_steps[0][i] + combined_steps[1][j];
                    results.push_back(lhs + combined_steps[2][k]);
                }


    return results;
}

std::vector<ConstrainedBivariatePolynomial<2>>
bottom_to_right_axis_paths(const Cell& cell, Line& axis) {
    
    auto results = std::vector<ConstrainedBivariatePolynomial<2>>();

    double sx = cell.s.x; //- cell.mid.x;
    double sy = cell.s.y; //- cell.mid.y;
    double tx = cell.t.x; //- cell.mid.x;
    double ty = cell.t.y; //- cell.mid.y;

    double alpha = angle(cell.s1, cell.t1);
    double beta = angle(cell.s2, cell.t2);

    Point s1 = cell.s1;
    Point t1 = cell.t1;
    Point s2 = cell.s2;
    Point t2 = cell.t2;

    double m = axis.grad_int().x;
    double b = axis.grad_int().y;

    auto h_to_axis_terms = horizontal_int(
        BivariatePolynomial<1>({{{{0, 0}},{{1, 0}}}}),
        BivariatePolynomial<1>({{{{sy/m-b/m, 0}},{{0, 0}}}}),
        cell, sy, true
    );

    auto v_to_axis_terms = vertical_int(
        BivariatePolynomial<1>({{{{sy, 0}},{{0, 0}}}}),
        BivariatePolynomial<1>({{{{b, 0}},{{m, 0}}}}),
        cell, 1, false
    );

    auto along_axis_from_x = axis_int(
        BivariatePolynomial<1>({{{{0, 0}},{{1, 0}}}}),
        BivariatePolynomial<1>({{{{-b/m, 0}},{{1/m, 0}}}}),
        cell, axis,
        {},
        {},
        {}

    );

    // auto alongaxis_from_int = axis_int(

    // );

    // auto h_from_axis_terms = horizontal_int(

    // );

    // auto v_from_axis_terms = vertical_int(

    // ); 
}

std::vector<ConstrainedBivariatePolynomial<2>>
bottom_right_axis_integrals(Line axis, const Cell& cell) {
    double sx = cell.s.x; //- cell.mid.x;
    double sy = cell.s.y; //- cell.mid.y;
    double tx = cell.t.x; //- cell.mid.x;
    double ty = cell.t.y; //- cell.mid.y;

    double m = axis.grad_int().x;
    double b = axis.grad_int().y;

    auto costs = std::vector<ConstrainedBivariatePolynomial<2>>();

    // up, axis, across

    // vertical(
    //  hi_lim, low_lim, cell,
    //  x_coeff, is_x_const
    // )

    auto uaa_up_terms = vertical_int(
        BivariatePolynomial<1>({{{{sy, 0}},{{0, 0}}}}),
        BivariatePolynomial<1>({{{{b, 0}},{{m, 0}}}}),
        cell, 1, false
    );

    double x_int = axis.getX(sy);
    double y_int = axis.getY(tx);

    // X, Y

    // lo = X
    // hi = 
    // Y = m*hi+b
    // y=mx+b
    // Y = mx+b

    auto uaa_axis_terms = axis_int(
        BivariatePolynomial<1>({{{{0, 0}},{{1, 0}}}}),
        BivariatePolynomial<1>({{{{-b/m, 1/m}},{{0, 0}}}}),
        cell, axis,
        {Polynomial<1>({sx, 0})}, 
        {Polynomial<1>({tx, 0}), Polynomial<1>({-b/m, 1/m})},
        {sy, sy ? std::min(ty, y_int) < sy : std::min(ty, y_int)}
    );

    auto uaa_across_terms = horizontal_int(
        BivariatePolynomial<1>({{{{-b/m, 1/m}},{{0, 0}}}}),
        BivariatePolynomial<1>({{{{tx, 0}},{{0, 0}}}}),
        cell, 1, false
    );


    auto up_axis_across = combine_steps(3, {uaa_up_terms, uaa_axis_terms, uaa_across_terms});
    for (auto poly: up_axis_across)
        costs.push_back(poly);
    // across, axis, up

    auto aau_across_terms = horizontal_int(
        BivariatePolynomial<1>({{{{sx, 0}},{{0, 0}}}}),
        BivariatePolynomial<1>({{{{(sy-b)/m, 0}},{{0, 0}}}}),
        cell, 1, false
    );

    auto aau_axis_terms = axis_int(
        BivariatePolynomial<1>({{{{(sy-b)/m, 0}},{{0, 0}}}}),
        BivariatePolynomial<1>({{{{tx, 0}},{{0, 0}}}}),
        cell, axis,
        {Polynomial<1>({sx, 0})}, 
        {Polynomial<1>({tx, 0}), Polynomial<1>({(sy-b)/m, 0})},
        {std::max(sy, m*tx+b), ty}
    );

    auto aau_up_terms = vertical_int(
        BivariatePolynomial<1>({{{{m*tx+b, 0}},{{0, 0}}}}),
        BivariatePolynomial<1>({{{{0, 1}},{{0, 0}}}}),
        cell, tx, true);

    auto across_axis_up = combine_steps(3, {aau_up_terms, aau_axis_terms, aau_across_terms});
    for (auto poly: across_axis_up)
        costs.push_back(poly);
    // up axis up

    auto uau_up_terms_2 = vertical_int(
        BivariatePolynomial<1>({{{{sy, 0}},{{0, 0}}}}),
        BivariatePolynomial<1>({{{{b, 0}},{{1, 0}}}}),
        cell, 1, false
    );

    auto uau_axis_terms = axis_int(
        BivariatePolynomial<1>({{{{(sy-b)/m, 0}},{{0, 0}}}}),
        BivariatePolynomial<1>({{{{tx, 0}},{{0, 0}}}}),
        cell, axis,
        {Polynomial<1>({sx, 0})}, 
        {Polynomial<1>({tx, 0}), Polynomial<1>({(sy-b)/m, 0})},
        {std::max(sy, m*tx+b), ty}
    );

    auto uau_up_terms_1 = vertical_int(
        BivariatePolynomial<1>({{{{m*tx+b, 0}},{{0, 0}}}}),
        BivariatePolynomial<1>({{{{0, 0}},{{1, 0}}}}),
        cell, 1, false
    );

    auto up_axis_up = combine_steps(3, {uau_up_terms_1, uau_axis_terms, uau_up_terms_2});
    for (auto poly: up_axis_up)
        costs.push_back(poly);

    // across axis across

    auto aaa_across_terms_1 = horizontal_int(
        BivariatePolynomial<1>({{{{sx, 0}},{{0, 0}}}}),
        BivariatePolynomial<1>({{{{(sy-b)/m, 0}},{{0, 0}}}}),
        cell, 1, false
    );

    auto aaa_axis_terms = axis_int(
        BivariatePolynomial<1>({{{{(sy-b)/m, 0}},{{0, 0}}}}),
        BivariatePolynomial<1>({{{{-b/m, 1/m}},{{0, 0}}}}),
        cell, axis,
        {Polynomial<1>({sx, 0})}, 
        {Polynomial<1>({tx, 0}), Polynomial<1>({-b/m, 1/m})},
        {sy, sy ? std::min(ty, y_int) < sy : std::min(ty, y_int)}
    );

    auto aaa_across_terms_2 = horizontal_int(
        BivariatePolynomial<1>({{{{-b/m, 1/m}},{{0, 0}}}}),
        BivariatePolynomial<1>({{{{tx, 0}},{{0, 0}}}}),
        cell, 1, false
    );

    auto across_axis_across = combine_steps(3, {aaa_across_terms_1, aaa_axis_terms, aaa_across_terms_2});
    for (auto poly: across_axis_across)
        costs.push_back(poly);

    for (auto cost: costs) {
        if (valid_point(cost, {tx, 0}) && !approx_zero(cost.f({tx, 0})))
            assert(false);
        else if (valid_point(cost, {tx, 0}))
            std::cout << "cost: " << cost.f({tx, 0}) << std::endl;
    }

    return costs;
}

std::vector<ConstrainedBivariatePolynomial<2>>
bottom_to_right_costs_2D(const Cell& cell) {
    // Coordinates of the cell in a system where cell.mid lies on (0, 0).
    // This makes it much easier to formulate the bivariate polynomials and constraints, but does require us to
    // translate the cell back such that (s.x, s.y) = (0, 0) using .translate_xy(-sx, -sy).
    double sx = cell.s.x; //- cell.mid.x;
    double sy = cell.s.y; //- cell.mid.y;
    double tx = cell.t.x; //- cell.mid.x;
    double ty = cell.t.y; //- cell.mid.y;

    // std::cout << "b2r\n";

    std::vector<ConstrainedBivariatePolynomial<2>> costs = std::vector<ConstrainedBivariatePolynomial<2>>();

    auto axes = get_axes(cell);

    double alpha = angle(cell.s1, cell.t1);
    double beta = angle(cell.s2, cell.t2);

    Point s1 = cell.s1;
    Point t1 = cell.t1;
    Point s2 = cell.s2;
    Point t2 = cell.t2;

    // solve_integral(BivariatePolynomial<1> low_lim, 
    // BivariatePolynomial<1> hi_lim,
    // BivariatePolynomial<1> integrand, double t_coeff,
    // Interval y_range, std::vector<Polynomial<1>> left_constraints, std::vector<Polynomial<1>> right_constraints)

    // if (approx_equal(cell.len1, 0.13726828480384534) && approx_equal(cell.len2, 0.18449462693244933)) {
    //     std::cout << "hi\n";
    // }
        // up and across

        auto ua_up_terms = vertical_int(
            BivariatePolynomial<1>({{{{sy, 0}},{{0, 0}}}}),
            BivariatePolynomial<1>({{{{0, 1}},{{0, 0}}}}),
            cell, 1, false
        );

        auto ua_across_terms = horizontal_int(
            BivariatePolynomial<1>({{{{0, 0}},{{1, 0}}}}),
            BivariatePolynomial<1>({{{{tx, 0}},{{0, 0}}}}),
            cell, 1, false
        );

        auto up_across = combine_steps(2, {ua_up_terms, ua_across_terms});

        for (auto poly: up_across)
            costs.push_back(poly);

        // across and up

        auto au_across_terms = horizontal_int(
            BivariatePolynomial<1>({{{{0, 0}},{{1, 0}}}}),
            BivariatePolynomial<1>({{{{tx, 0}},{{0, 0}}}}),
            cell, sy, true
        );

        auto au_up_terms = vertical_int(
            BivariatePolynomial<1>({{{{sy, 0}},{{0, 0}}}}),
            BivariatePolynomial<1>({{{{0, 1}},{{0, 0}}}}),
            cell, tx, true
        );

        auto across_up = combine_steps(2, {au_across_terms, au_up_terms});

        for (auto poly: across_up)
            costs.push_back(poly);

    if (axes.size() >= 1) {
        Line axis = axes[0];
        auto axis_costs = bottom_right_axis_integrals(axis, cell);
        for (auto poly: axis_costs)
            costs.push_back(poly);
    }

    if (axes.size() == 2) {
        Line axis = axes[1];
        auto axis_costs = bottom_right_axis_integrals(axis, cell);
        for (auto poly: axis_costs)
            costs.push_back(poly);
    }

    auto sus_polynomials = std::vector<ConstrainedBivariatePolynomial<2>>();

    
    for (auto cost: costs)
        if (approx_equal(cost.y_interval.max, ty))
            sus_polynomials.push_back(cost);

    return costs;
}

std::vector<ConstrainedBivariatePolynomial<2>>
vertical_int_b2t(BivariatePolynomial<1> low_lim, BivariatePolynomial<1> hi_lim,
const Cell& cell, std::string fixed_var) {

    auto sx = 0.;
    auto tx = cell.len1;
    auto sy = 0.;
    auto ty = cell.len1;

    double alpha = angle(cell.s1, cell.t1);
    double beta = angle(cell.s2, cell.t2);

    Point s1 = cell.s1;
    Point t1 = cell.t1;
    Point s2 = cell.s2;
    Point t2 = cell.t2;

    auto results = std::vector<ConstrainedBivariatePolynomial<2>>();
    std::vector<ConstrainedBivariatePolynomial<2>> x_terms;
    std::vector<ConstrainedBivariatePolynomial<2>> y_terms;

    if (fixed_var == "y") {

        x_terms = solve_integral(low_lim, hi_lim,
            BivariatePolynomial<1>({{{{s1.x-s2.x, cos(alpha)}},{{0, 0}}}}),
            cos(beta), {sy, ty},
            {Polynomial<1>({sx, 0})}, {Polynomial<1>({tx, 0})} 
        );

        y_terms = solve_integral(low_lim, hi_lim,
            BivariatePolynomial<1>({{{{s1.y-s2.y, sin(alpha)}},{{0, 0}}}}),
            sin(beta), {sy, ty},
            {Polynomial<1>({sx, 0})}, {Polynomial<1>({tx, 0})} 
        );

    } else if (fixed_var == "x") {

        x_terms = solve_integral(low_lim, hi_lim,
            BivariatePolynomial<1>({{{{s1.x-s2.x, 0}},{{cos(alpha), 0}}}}),
            cos(beta), {sy, ty},
            {Polynomial<1>({sx, 0})}, {Polynomial<1>({tx, 0})} 
        );

        y_terms = solve_integral(low_lim, hi_lim,
            BivariatePolynomial<1>({{{{s1.y-s2.y, 0}},{{sin(alpha), 0}}}}),
            sin(beta), {sy, ty},
            {Polynomial<1>({sx, 0})}, {Polynomial<1>({tx, 0})} 
        );

    }


    for (auto xterm: x_terms)
        for (auto yterm: y_terms)
            results.push_back(xterm + yterm);

    auto valid_results = std::vector<ConstrainedBivariatePolynomial<2>>();

    for (auto f: results)
        if (valid_constraints(f)) {
            if (!is_positive(f))
                std::cout << "nope\n";
            valid_results.push_back(f);
        }

    return valid_results;
}

std::vector<ConstrainedBivariatePolynomial<2>>
horizontal_int_b2t(BivariatePolynomial<1> low_lim, BivariatePolynomial<1> hi_lim,
const Cell& cell, double height) {

    auto sx = 0.;
    auto tx = cell.len1;
    auto sy = 0.;
    auto ty = cell.len1;

    double alpha = angle(cell.s1, cell.t1);
    double beta = angle(cell.s2, cell.t2);

    Point s1 = cell.s1;
    Point t1 = cell.t1;
    Point s2 = cell.s2;
    Point t2 = cell.t2;

    auto results = std::vector<ConstrainedBivariatePolynomial<2>>();

    auto x_terms = solve_integral(
        low_lim, hi_lim,
        BivariatePolynomial<1>({{{{s2.x-s1.x + height*cos(beta), 0}},{{0, 0}}}}),
        cos(alpha), {sy, ty},
        {Polynomial<1>({sx, 0})}, {Polynomial<1>({tx, 0})}
    );

    auto y_terms = solve_integral(
        low_lim, hi_lim,
        BivariatePolynomial<1>({{{{s2.y-s1.y + height*sin(beta), 0}},{{0, 0}}}}),
        sin(alpha), {sy, ty},
        {Polynomial<1>({sx, 0})}, {Polynomial<1>({tx, 0})}
    );

    for (auto xterm: x_terms)
        for (auto yterm: y_terms)
            results.push_back(xterm + yterm);

    auto valid_results = std::vector<ConstrainedBivariatePolynomial<2>>();

    for (auto f: results)
        if (valid_constraints(f)) {
            if (!is_positive(f))
                std::cout << "nope\n";
            valid_results.push_back(f);
        }

    return valid_results;
}

std::vector<ConstrainedBivariatePolynomial<2>>
axis_int_b2t(BivariatePolynomial<1> low_lim, BivariatePolynomial<1> hi_lim, const Cell& cell, Line axis,
std::vector<Polynomial<1>> left_constraints, std::vector<Polynomial<1>> right_constraints, 
Interval y_range) {
    Point s1 = cell.s1;
    Point s2 = cell.s2;
    Point t1 = cell.t1;
    Point t2 = cell.t2;

    double sx = cell.s.x; //- cell.mid.x;
    double sy = cell.s.x; //- cell.mid.y;
    double tx = cell.t.x; //- cell.mid.x;
    double ty = cell.t.x; //- cell.mid.y;

    double alpha = angle(s1, t1);
    double beta = angle(s2, t2);

    auto gi = axis.grad_int();

    double m = gi.x;
    double b = gi.y;

    auto xterms = solve_integral(
        low_lim, hi_lim,
        BivariatePolynomial<1>({{{{s1.x-s2.x - b*cos(beta), 0}},{{0, 0}}}}),
        m*cos(beta) - cos(alpha), y_range,
        left_constraints, right_constraints
    );

    auto yterms = solve_integral(
        low_lim, hi_lim,
        BivariatePolynomial<1>({{{{s1.y-s2.y - b*sin(beta), 0}},{{0, 0}}}}),
        m*sin(beta) - sin(alpha), y_range,
        left_constraints, right_constraints
    );

    std::vector<ConstrainedBivariatePolynomial<2>> results = std::vector<ConstrainedBivariatePolynomial<2>>();

    for (auto xterm: xterms)
        for (auto yterm: yterms)
            results.push_back(xterm + yterm);

    auto valid_results = std::vector<ConstrainedBivariatePolynomial<2>>();

    for (auto f: results)
        if (valid_constraints(f))
            valid_results.push_back(f);

    return valid_results;
}

std::vector<ConstrainedBivariatePolynomial<2>>
bottom_top_axis_integrals(Line axis, const Cell& cell) {

    auto sx = 0.;
    auto tx = cell.len1;
    auto sy = 0.;
    auto ty = cell.len1;

    double alpha = angle(cell.s1, cell.t1);
    double beta = angle(cell.s2, cell.t2);

    Point s1 = cell.s1;
    Point t1 = cell.t1;
    Point s2 = cell.s2;
    Point t2 = cell.t2;

    double m = axis.grad_int().x;
    double b = axis.grad_int().y;

    auto results = std::vector<ConstrainedBivariatePolynomial<2>>();

    auto bottom_int = axis.getX(0);
    auto top_int = axis.getX(cell.len2);

    // up-axis-up

    auto uau_u1_terms = vertical_int_b2t(
        BivariatePolynomial<1>({{{{0, 0}},{{0, 0}}}}),
        BivariatePolynomial<1>({{{{b, 0}},{{m, 0}}}}),
        cell, "x"
    );

    auto uau_axis_terms = axis_int(
        BivariatePolynomial<1>({{{{0, 0}},{{1, 0}}}}),
        BivariatePolynomial<1>({{{{0, 1}},{{0, 0}}}}),
        cell, axis,
        {Polynomial<1>({sx, 0}), Polynomial<1>({bottom_int})}, 
        {Polynomial<1>({tx, 0})},
        {sy, std::min(ty, top_int)}
    );

    auto uau_u2_terms = vertical_int_b2t(
        BivariatePolynomial<1>({{{{b, 0}},{{m, 0}}}}),
        BivariatePolynomial<1>({{{{cell.len2 , 0}},{{0, 0}}}}),
        cell, "y"
    );

    auto up_axis_up = combine_steps(3, {uau_u1_terms, uau_axis_terms, uau_u2_terms});
    for (auto poly: up_axis_up)
        results.push_back(poly);


    // across-axis-up

    auto aau_across_terms = horizontal_int_b2t(
        BivariatePolynomial<1>({{{{0, 0}},{{1, 0}}}}),
        BivariatePolynomial<1>({{{{bottom_int, 0}},{{0, 0}}}}),
        cell, 0
    );

    auto aau_axis_terms = axis_int(
        BivariatePolynomial<1>({{{{bottom_int, 0}},{{0, 0}}}}),
        BivariatePolynomial<1>({{{{0, 1}},{{0, 0}}}}),
        cell, axis,
        {Polynomial<1>({sx, 0})}, 
        {Polynomial<1>({tx, 0}), Polynomial<1>({bottom_int, 0})},
        {sy, std::min(ty, top_int)}
    );

    auto aau_up_terms = vertical_int_b2t(
        BivariatePolynomial<1>({{{{b, 0}},{{m, 0}}}}),
        BivariatePolynomial<1>({{{{cell.len2 , 0}},{{0, 0}}}}),
        cell, "y"
    );

    auto across_axis_up = combine_steps(3, {aau_across_terms, aau_axis_terms, aau_up_terms});
    for (auto poly: across_axis_up)
        results.push_back(poly);

    // up-axis-across

    auto uaa_up_terms = vertical_int_b2t(
        BivariatePolynomial<1>({{{{0, 0}},{{0, 0}}}}),
        BivariatePolynomial<1>({{{{b, 0}},{{m, 0}}}}),
        cell, "x"
    );

    auto uaa_axis_terms = axis_int(
        BivariatePolynomial<1>({{{{0, 0}},{{1, 0}}}}),
        BivariatePolynomial<1>({{{{top_int, 0}},{{0, 0}}}}),
        cell, axis,
        {Polynomial<1>({sx, 0})}, 
        {Polynomial<1>({tx, 0}), Polynomial<1>({bottom_int, 0})},
        {std::max(sy, top_int), ty}
    );

    auto uaa_across_terms = horizontal_int_b2t(
        BivariatePolynomial<1>({{{{top_int, 0}},{{0, 0}}}}),
        BivariatePolynomial<1>({{{{0, 1}},{{0, 0}}}}),
        cell, cell.len2
    );

    auto up_axis_across = combine_steps(3, {uaa_up_terms, uaa_axis_terms, uaa_across_terms});
    for (auto poly: across_axis_up)
        results.push_back(poly);

    // across-axis-across

    auto aaa_a1_terms = horizontal_int_b2t(
        BivariatePolynomial<1>({{{{0, 0}},{{1, 0}}}}),
        BivariatePolynomial<1>({{{{bottom_int, 0}},{{0, 0}}}}),
        cell, 0
    );

    auto aaa_axis_terms = axis_int(
        BivariatePolynomial<1>({{{{bottom_int, 0}},{{0, 0}}}}),
        BivariatePolynomial<1>({{{{top_int, 0}},{{0, 0}}}}),
        cell, axis,
        {Polynomial<1>({sx, 0})}, 
        {Polynomial<1>({tx, 0}), Polynomial<1>({bottom_int, 0})},
        {std::max(sy, top_int), ty}
    );

    auto aaa_a2_terms = horizontal_int_b2t(
        BivariatePolynomial<1>({{{{top_int, 0}},{{0, 0}}}}),
        BivariatePolynomial<1>({{{{0, 1}},{{0, 0}}}}),
        cell, cell.len2
    );

    auto across_axis_across = combine_steps(3, {aaa_a1_terms, aaa_axis_terms, aaa_a2_terms});
    for (auto poly: across_axis_up)
        results.push_back(poly);

    auto valid_results = std::vector<ConstrainedBivariatePolynomial<2>>();

    for (auto f: results)
        if (valid_constraints(f))
            valid_results.push_back(f);

    return valid_results;
}


std::vector<ConstrainedBivariatePolynomial<2>>
bottom_to_top_costs_2D(const Cell& cell) {
    double sx = cell.s.x; //- cell.mid.x;
    double sy = cell.s.x; //- cell.mid.y;
    double tx = cell.t.x; //- cell.mid.x;
    double ty = cell.t.x; //- cell.mid.y;

    // std::cout << "b2t\n";

    std::vector<ConstrainedBivariatePolynomial<2>> costs = std::vector<ConstrainedBivariatePolynomial<2>>();

    auto axes = get_axes(cell);

    double alpha = angle(cell.s1, cell.t1);
    double beta = angle(cell.s2, cell.t2);

    Point s1 = cell.s1;
    Point t1 = cell.t1;
    Point s2 = cell.s2;
    Point t2 = cell.t2;

    // up and across

    auto ua_up_terms = vertical_int_b2t(
        BivariatePolynomial<1>({{{{0, 0}},{{0, 0}}}}),
        BivariatePolynomial<1>({{{{cell.len2, 0}},{{0, 0}}}}),
        cell, "x"
    );
    
    auto ua_across_terms = horizontal_int_b2t(
        BivariatePolynomial<1>({{{{0, 0}},{{1, 0}}}}),
        BivariatePolynomial<1>({{{{0, 1}},{{0, 0}}}}),
        cell, cell.len2
    );

    auto up_across = combine_steps(2, {ua_across_terms, ua_up_terms});
    for (auto poly: up_across) {
        if (!is_positive(poly))
                std::cout << "nope\n";
        costs.push_back(poly);
    }

    //across and up

    auto au_across_terms = horizontal_int_b2t(
        BivariatePolynomial<1>({{{{0, 0}},{{1, 0}}}}),
        BivariatePolynomial<1>({{{{0, 1}},{{0, 0}}}}),
        cell, 0
    );

    auto au_up_terms = vertical_int_b2t(
        BivariatePolynomial<1>({{{{0, 0}},{{0, 0}}}}),
        BivariatePolynomial<1>({{{{cell.len2, 0}},{{0, 0}}}}),
        cell, "y"
    );

    auto across_up = combine_steps(2, {au_across_terms, au_up_terms});
    for (auto poly: across_up) {
        if (!is_positive(poly))
                std::cout << "nope\n";
        costs.push_back(poly);
    }


    if (axes.size() == 1) {
        Line axis = axes[0];
        auto axis_costs = bottom_top_axis_integrals(axis, cell);
        for (auto poly: axis_costs) {
            if (!is_positive(poly))
                std::cout << "nope\n";
            costs.push_back(poly);
        }
    }

    if (axes.size() == 2) {
        Line axis = axes[1];
        auto axis_costs = bottom_top_axis_integrals(axis, cell);
        for (auto poly: axis_costs) {
            if (!is_positive(poly))
                std::cout << "nope\n";
            costs.push_back(poly);
        }
    }

    // std::cout << "b2t cost count: " << costs.size() << std::endl;
    
    return costs;
}

template<class Iterator>
PiecewisePolynomial<2> propagate(
    const std::vector<ConstrainedBivariatePolynomial<2>>& cell_costs,
    Iterator pieces_it,
    Iterator pieces_end
) {
    std::vector<PolynomialPiece<2>> min_pieces;

        #ifndef NDEBUG
        std::vector<ConstrainedBivariatePolynomial<2>> total_costs;
        #endif
        int counter = 0;
        PiecewisePolynomial<2> out_cost;
        while (pieces_it != pieces_end) {
            // piece_in_cost: cost of optimal path from origin to point on in-boundary
            const PolynomialPiece<2>& piece_in_cost = *pieces_it;
            counter++;
            // piece_out_cost: cost of optimal path from origin to point on out-boundary
            PiecewisePolynomial<2> piece_out_cost;

            // cell_cost: cost of optimal path from point on in-boundary to point on out-boundary
            for (const auto& cell_cost : cell_costs) {
                // total_cost: cost of optimal path from origin through point on in-boundary to point on out-boundary
                const auto total_cost = cell_cost.add_x(piece_in_cost);

                #ifndef NDEBUG
                total_costs.push_back(total_cost);
                #endif

                const PiecewisePolynomial<2> min_total_cost = find_minimum(
                    total_cost.f,
                    total_cost.y_interval,
                    total_cost.left_constraints,
                    total_cost.right_constraints
                );

                min_pieces.insert(min_pieces.end(), min_total_cost.pieces.begin(), min_total_cost.pieces.end());
            }
            std::cout << "incoming pieces: " << counter << std::endl;
            ++pieces_it;
        }

        const auto result = naive_lower_envelope(min_pieces);

        #ifndef NDEBUG
        verify_minimum(result, total_costs);
        #endif

        return result;
}
}

/**
 * Compute the cost function from the origin of the given cell to a point A along the bottom boundary of the cell.
 * Note that you can compute a similar function along the left boundary by providing a transposed cell.
 *
 * Used as the base case in the CDTW dynamic program.
 *
 * @param cell
 * @return Piecewise polynomial over the domain 0 <= A <= cell.width.
 */
template<>
PiecewisePolynomial<2> CDTW<2, Norm::L1, Norm::L1>::base_bottom(const Cell& cell) const {
    // TODO: Clean this up; sx, sy, tx, ty are assumed to be coordinates in a space where (0,0) is the ellipse center
    double sx = cell.s.x; //- cell.mid.x;
    double sy = cell.s.y; //- cell.mid.y;
    double tx = cell.t.x; //- cell.mid.x;
    double ty = cell.t.y;

    auto costs = horizontal_int_b2t(
        BivariatePolynomial<1>({{{{0, 0}},{{0, 0}}}}),
        BivariatePolynomial<1>({{{{0, 0}},{{1, 0}}}}),
        cell, 0
    );

    auto results = PiecewisePolynomial<2>();

    for (auto cost: costs) {
        auto poly = cost.f.slice_at_y(0);
        double max = cell.len1;
        double min = 0;
        for (auto rc: cost.right_constraints) 
            if (rc(0) < max)
                max = rc(0);
        
        for (auto lc: cost.left_constraints)
            if (lc(0) > min)
                min = lc(0);

        if (max > min)
            results.pieces.push_back({{min, max}, poly});
    }
    
    auto final_result = naive_lower_envelope(results.pieces);

    return final_result;
    // bool same_direction = (cell.s1.x <= cell.t1.x) == (cell.s2.x <= cell.t2.x);

    // if (!same_direction) { // Downwards valley
    //     if (sx >= -sy) {
    //         return PiecewisePolynomial<2>({
    //             {{0, tx-sx}, Polynomial<2>({0, sx + sy, 1./2})},
    //         });
    //     } else if (-sy >= tx) {
    //         return PiecewisePolynomial<2>({
    //             {{0, tx-sx}, Polynomial<2>({0, -sx-sy, -1./2})},
    //         });
    //     } else { // sx < -sy < tx
    //         return PiecewisePolynomial<2>({
    //             {{0, -sy-sx}, Polynomial<2>({0, -sx-sy, -1./2})},
    //             {{-sy-sx, tx-sx}, Polynomial<2>({(sx+sy) * (sx+sy), sx+sy, 1./2})},
    //         });
    //     }
    // }

    // // Upwards valley

    // if (sy <= sx) {
    //     return PiecewisePolynomial<2>({
    //         {{sx, tx}, Polynomial<2>({-(sx * sx) / 2 + sx * sy, -sy, 1. / 2})}
    //     }).translate(-sx);
    // } else if (sy >= tx) {
    //     return PiecewisePolynomial<2>({
    //         {{sx, tx}, Polynomial<2>({(sx * sx) / 2 - sx * sy, sy, -1. / 2})}
    //     }).translate(-sx);
    // } else { // sx < sy < tx
    //     return PiecewisePolynomial<2>({
    //         {{sx, sy}, Polynomial<2>({(sx * sx) / 2 - sx * sy, sy, -1. / 2})},
    //         {{sy, tx}, Polynomial<2>({(sx * sx) / 2 - sx * sy + sy * sy, -sy, 1. / 2})}
    //     }).translate(-sx);
    // }
    
    
    return PiecewisePolynomial<2>();
}

/**
 * Compute the piecewise bivariate polynomial representing the cost of the optimal path between a point A on the bottom
 * boundary and a point B on the right boundary.
 *
 * @param cell
 * @return A piecewise bivariate polynomial over the domain 0 <= A <= cell.width and 0 <=B <= cell.height.
 */
template<>
std::vector<ConstrainedBivariatePolynomial<2>>
CDTW<2, Norm::L1, Norm::L1>::bottom_to_right_costs(const Cell& cell) const {
    // auto subdivision = get_subdivision(cell);

    // if (subdivision.size() == 0)
    //     return bottom_to_right_costs_2D(cell);

    // auto& bl = subdivision[0];
    // auto& tl = subdivision[1];
    // auto& br = subdivision[2];
    // auto& tr = subdivision[3];

    // // auto& br_in_functions

    // std::cout << subdivision.size() << std::endl;

    auto costs = bottom_to_right_costs_2D(cell);

    // for (auto cost: costs)
    //     if (!is_positive(cost))
    //         std::cout << "...\n";

    return costs;
}

/**
 * Compute the piecewise bivariate polynomial representing the cost of the optimal path between a point A on the bottom
 * boundary and a point B on the top boundary.
 *
 * @param cell
 * @return A piecewise bivariate polynomial over the domain 0 <= a <= b <= cell.width.
 */
template<>
std::vector<ConstrainedBivariatePolynomial<2>>
CDTW<2, Norm::L1, Norm::L1>::bottom_to_top_costs(const Cell& cell) const {
    auto subdivision = get_subdivision(cell);

    auto result = bottom_to_top_costs_2D(cell);

    for (auto cost: result)
        if (!is_positive(cost))
            std::cout << "...\n";

    if (subdivision.size() == 0)
        return bottom_to_top_costs_2D(cell);

    return bottom_to_top_costs_2D(cell);
}

#endif //TRAJECTORY_CLUSTERING_2D_L1_L1_H