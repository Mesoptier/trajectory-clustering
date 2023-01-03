#pragma once

#include "ConstrainedBivariatePolynomial.h"
#include "2d_utils.h"

std::vector<Polynomial<1>> x_constraints(
    BivariatePolynomial<1> lhs, BivariatePolynomial<1> rhs, 
    std::string ineq, std::string side);


std::vector<double> y_constraints(
    BivariatePolynomial<1> lhs, BivariatePolynomial<1> rhs, 
    std::string ineq, std::string side);

void update_constraints(BivariatePolynomial<1> lhs, BivariatePolynomial<1> rhs, std::string ineq,
    std::vector<Polynomial<1>>& left_constraints, std::vector<Polynomial<1>>& right_constraints,
    double& y_lower_bound, double& y_upper_bound);

/**
 * Explicity solves an integral and returns a vector
 * of constrained bivariate quadratic functions
*/
std::vector<ConstrainedBivariatePolynomial<2>>
solve_integral(BivariatePolynomial<1> low_lim, 
BivariatePolynomial<1> hi_lim,
BivariatePolynomial<1> integrand, double t_coeff,
Interval_c y_range, std::vector<Polynomial<1>> left_constraints, std::vector<Polynomial<1>> right_constraints);