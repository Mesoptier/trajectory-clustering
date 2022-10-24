#ifndef TRAJECTORY_CLUSTERING_2D_L1_L1_H
#define TRAJECTORY_CLUSTERING_2D_L1_L1_H

#include "cdtw.h"
#include "BivariatePolynomial.h"
// #include "1d-l1-l1.h"

//
// 2D + L1 image norm + L1 param norm
//

// struct Constraint {
//     bool is_y_constraint;
//     bool is_lower_bound;
//     double y_bound;
//     Polynomial<1> x_bound;

//     Constraint(bool y_const, bool lower_bound, double _y_bound, Polynomial<1> _x_bound) :
//     is_y_constraint(y_const), is_lower_bound(lower_bound), y_bound(_y_bound), x_bound(_x_bound) {}

//     bool valid_point(Point p);

//     Point normal();

//     void shift();

//     Point boundary_point();

//     bool zero_between(const Constraint& other);

//     bool operator==(const Constraint& rhs);

//     bool operator!=(const Constraint& rhs);
// };

// bool domain_covered(std::vector<ConstrainedBivariatePolynomial<2>> costs, const Cell& cell);

// bool is_parallel(Constraint a, Constraint b);

// Point intersect(Constraint a, Constraint b);

// bool valid_triple(Constraint a, Constraint b, Constraint c);

// bool valid_constraints(ConstrainedBivariatePolynomial<2>& poly);

/**
 * 
 * Samples points from the domian
 * of a constrained bivariate quadratic
 */
std::vector<Point> sample_domain(ConstrainedBivariatePolynomial<2>& poly);


bool valid_point(ConstrainedBivariatePolynomial<2>& poly, Point p);

bool is_positive(ConstrainedBivariatePolynomial<2>& poly);



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
Interval_c y_range, std::vector<Polynomial<1>> left_constraints, std::vector<Polynomial<1>> right_constraints);

//TODO: Implement this properly
bool has_positive_slope(Line line);

distance_t get_angle(Line line);

bool contains_midpoint(const Cell& cell);


std::vector<Cell> 
get_subdivision(const Cell& cell);

double param_space_coord(const Point& s, const Point& t, Point a);

Point param_space(const Point& s1, const Point& t1,
const Point& s2, const Point& t2,
Point a, Point b);

/**
 * Determines if axis is monote
 * @param axis
 * @return true if axis is monotone, false otherwise
 * */
bool is_monotone(Line& axis);

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
CellIntersections get_intersections(const Cell& cell, Line& line);

bool intersects_cell(const Cell& cell, CellIntersections& intersections);

Line axis_from_points(const Cell& cell, 
Point l1_p1, Point l2_p1, Point l1_p2, Point l2_p2);


/**
* Get axes for given cell
* @param cell
* @return vector of lines corresponding to the axes which
* intersect the cell
*/
std::vector<Line>
get_axes(const Cell& cell);


std::vector<ConstrainedBivariatePolynomial<2>>
horizontal_int(BivariatePolynomial<1> low_lim, BivariatePolynomial<1> hi_lim,
const Cell& cell, double y_coef, bool y_const);

std::vector<ConstrainedBivariatePolynomial<2>>
vertical_int(BivariatePolynomial<1> low_lim, BivariatePolynomial<1> hi_lim,
const Cell& cell, double x_coef, bool x_const
);


std::vector<ConstrainedBivariatePolynomial<2>>
axis_int(BivariatePolynomial<1> low_lim, BivariatePolynomial<1> hi_lim, const Cell& cell, Line axis,
std::vector<Polynomial<1>> left_constraints, std::vector<Polynomial<1>> right_constraints, 
Interval_c y_range);

std::vector<ConstrainedBivariatePolynomial<2>>
combine_steps(size_t step_count, std::vector<std::vector<ConstrainedBivariatePolynomial<2>>> steps);

std::vector<ConstrainedBivariatePolynomial<2>>
combine_terms(size_t step_count, 
std::vector<std::vector<ConstrainedBivariatePolynomial<2>>> xterms,
std::vector<std::vector<ConstrainedBivariatePolynomial<2>>> yterms);

std::vector<ConstrainedBivariatePolynomial<2>>
bottom_to_right_axis_paths(const Cell& cell, Line& axis);

std::vector<ConstrainedBivariatePolynomial<2>>
bottom_right_axis_integrals(Line axis, const Cell& cell);

std::vector<ConstrainedBivariatePolynomial<2>>
bottom_to_right_costs_2D(const Cell& cell);

std::vector<ConstrainedBivariatePolynomial<2>>
vertical_int_b2t(BivariatePolynomial<1> low_lim, BivariatePolynomial<1> hi_lim,
const Cell& cell, std::string fixed_var);

std::vector<ConstrainedBivariatePolynomial<2>>
horizontal_int_b2t(BivariatePolynomial<1> low_lim, BivariatePolynomial<1> hi_lim,
const Cell& cell, double height);

std::vector<ConstrainedBivariatePolynomial<2>>
axis_int_b2t(BivariatePolynomial<1> low_lim, BivariatePolynomial<1> hi_lim, const Cell& cell, Line axis,
std::vector<Polynomial<1>> left_constraints, std::vector<Polynomial<1>> right_constraints, 
Interval_c y_range);


std::vector<ConstrainedBivariatePolynomial<2>>
bottom_top_axis_integrals(Line axis, const Cell& cell);

std::vector<ConstrainedBivariatePolynomial<2>>
bottom_to_top_costs_2D(const Cell& cell);

template<class Iterator>
PiecewisePolynomial<2> propagate(
    const std::vector<ConstrainedBivariatePolynomial<2>>& cell_costs,
    Iterator pieces_it,
    Iterator pieces_end
);


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
PiecewisePolynomial<2> CDTW<2, Norm::L1, Norm::L1>::base_bottom(const Cell& cell) const;

/**
 * Compute the piecewise bivariate polynomial representing the cost of the optimal path between a point A on the bottom
 * boundary and a point B on the right boundary.
 *
 * @param cell
 * @return A piecewise bivariate polynomial over the domain 0 <= A <= cell.width and 0 <=B <= cell.height.
 */
template<>
std::vector<ConstrainedBivariatePolynomial<2>>
CDTW<2, Norm::L1, Norm::L1>::bottom_to_right_costs(const Cell& cell) const;
/**
 * Compute the piecewise bivariate polynomial representing the cost of the optimal path between a point A on the bottom
 * boundary and a point B on the top boundary.
 *
 * @param cell
 * @return A piecewise bivariate polynomial over the domain 0 <= a <= b <= cell.width.
 */
template<>
std::vector<ConstrainedBivariatePolynomial<2>>
CDTW<2, Norm::L1, Norm::L1>::bottom_to_top_costs(const Cell& cell) const;

#endif //TRAJECTORY_CLUSTERING_2D_L1_L1_H