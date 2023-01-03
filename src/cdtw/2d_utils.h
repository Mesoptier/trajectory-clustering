#pragma once
#include "../IntegralFrechet/Cell.h"
#include "ConstrainedBivariatePolynomial.h"

namespace utils {

long double angle(Point s, Point t);

double _sin(double angle);

double _cos(double angle);

enum VARIABLE_TYPE { X, Y, _ };

struct Constraint {
    bool is_y_constraint;
    bool is_lower_bound;
    double y_bound;
    Polynomial<1> x_bound;

    Constraint(bool y_const, bool lower_bound, double _y_bound, Polynomial<1> _x_bound) :
    is_y_constraint(y_const), is_lower_bound(lower_bound), y_bound(_y_bound), x_bound(_x_bound) {}

    /**
     * @brief check that p satisfies the constraint
     * 
     * @param p 
     * @return true 
     * @return false 
     */
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


    /**
     * @brief returns a point representing a vector normal
     * to the constraint
     * 
     * @return Point 
     */
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

    bool zero_between(const Constraint& other) {
        if (!is_y_constraint && !other.is_y_constraint) 
            if (is_lower_bound && !other.is_lower_bound || !is_lower_bound && other.is_lower_bound) 
                if (x_bound == other.x_bound)
                    return true;

        if (is_y_constraint && other.is_y_constraint) 
            if (is_lower_bound && !other.is_lower_bound || !is_lower_bound && other.is_lower_bound) 
                if (approx_equal(y_bound, other.y_bound))
                    return true;
            
        return false;
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

/**
 * @brief Check that every point in the cell is within the
 * region of at least one polynomial
 * 
 * @param costs 
 * @param cell 
 * @return true 
 * @return false 
 */
bool domain_covered(std::vector<ConstrainedBivariatePolynomial<2>> costs, const Cell& cell);


bool is_parallel(Constraint a, Constraint b);


Point intersect(Constraint a, Constraint b);
/**
 * @brief Check that a, b and c define a non-empty region.
 * 
 * @param a 
 * @param b 
 * @param c 
 * @return true 
 * @return false 
 */
bool valid_triple(Constraint a, Constraint b, Constraint c);
bool valid_constraints(ConstrainedBivariatePolynomial<2>& poly);

/**
 * Samples points from the domian
 * of a constrained bivariate quadratic
 */
std::vector<Point> sample_domain(ConstrainedBivariatePolynomial<2>& poly);

bool valid_point(ConstrainedBivariatePolynomial<2>& poly, Point p);

bool is_positive(ConstrainedBivariatePolynomial<2>& poly);

bool contains_midpoint(const Cell& cell);

std::vector<Cell> 
get_subdivision(const Cell& cell);

double param_space_coord(const Point& s, const Point& t, Point a);

Point param_space(const Point& s1, const Point& t1,
const Point& s2, const Point& t2,
Point a, Point b);

bool is_monotone(Line& axis);

bool intersect(const Cell& cell, Line v);

Line axis_from_point_dir(const Cell& cell, Point l1_p, Point l2_p, Point dir);

Line axis_from_points(const Cell& cell, 
Point l1_p1, Point l2_p1, Point l1_p2, Point l2_p2);


distance_t l1_h_function(const Cell& cell, Point p);

/*
* Returns the proper axis of the arugument cell.
* v1 and v2 must be the axes of the cell.
*/
Line get_proper_axis(const Cell& cell, Line v1, Line v2);

/**
* Get axes for given cell
*/
std::vector<Line>
get_axis(const Cell& cell);

}