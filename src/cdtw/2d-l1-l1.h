#ifndef TRAJECTORY_CLUSTERING_2D_L1_L1_H
#define TRAJECTORY_CLUSTERING_2D_L1_L1_H

#include "cdtw.h"
#include "BivariatePolynomial.h"


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
bottom_right_valley_integrals(Line axis, const Cell& cell);

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
bottom_top_valley_integrals(Line axis, const Cell& cell);

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