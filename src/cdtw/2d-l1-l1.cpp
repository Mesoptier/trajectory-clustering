#ifndef TRAJECTORY_CLUSTERING_2D_L1_L1_H
#define TRAJECTORY_CLUSTERING_2D_L1_L1_H

#include "cdtw.h"
#include "BivariatePolynomial.h"
#include "2d-l1-l1.h"
#include "2d_utils.h"
#include "integral_solver.h"

//
// 2D + L1 image norm + L1 param norm
//

using namespace utils;

namespace {


using Cnstrnts = std::vector<Polynomial<1>>;

/**
 * @brief Computes a vector of constrained bivariate polynomials
 * representing the cost of integrating along a horizontal line
 * from low_lim to hi_lim.
 * 
 * @param low_lim 
 * @param hi_lim 
 * @param cell 
 * @param y_coef 
 * @param y_const Determines whether y is to be treated as a constant.
 * @return std::vector<ConstrainedBivariatePolynomial<2>> 
 */
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

    double sx = cell.s.x;
    double sy = cell.s.y;
    double tx = cell.t.x;
    double ty = cell.t.y; 

    double alpha = angle(s1, t1);
    double beta = angle(s2, t2);

    if (y_const) {
        auto terms_x = solve_integral(low_lim, hi_lim,
            // Polynomial
            BivariatePolynomial<1>({{{{s1.x-s2.x - y_coef*_cos(beta), 0}},{{0, 0}}}}),
            // a
            -_cos(alpha), 
            {sy, ty},
            {Polynomial<1>({sx, 0})}, {Polynomial<1>({tx, 0})} 
        );

        for (auto term: terms_x)  
            xterms.push_back(term);

        auto terms_y = solve_integral(low_lim, hi_lim,
            // polynomial
            BivariatePolynomial<1>({{{{s1.y-s2.y - y_coef*_sin(beta), 0}},{{0, 0}}}}),
            // a
            -_sin(alpha), 
            {sy, ty},
            {Polynomial<1>({sx, 0})}, {Polynomial<1>({tx, 0})} 
        );

        for (auto term: terms_y)  
            yterms.push_back(term);

            
    } else {
        auto terms_x = solve_integral(low_lim, hi_lim,
            // polynomial
            BivariatePolynomial<1>({{{{s1.x-s2.x, -_cos(beta)}},{{0, 0}}}}),
            // a
            -_cos(alpha), 
            {sy, ty},
            {Polynomial<1>({sx, 0})}, {Polynomial<1>({tx, 0})} 
        );

        for (auto term: terms_x)  
            xterms.push_back(term);

        auto terms_y = solve_integral(low_lim, hi_lim,
            // polynomial
            BivariatePolynomial<1>({{{{s1.y-s2.y, -_sin(beta)}},{{0, 0}}}}),
            // a
            -_sin(alpha), 
            {sy, ty},
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
        if (valid_constraints(f)) {
            valid_results.push_back(f);
        }

    return valid_results;
}

/**
 * @brief Computes a vector of constrained bivariate polynomials
 * representing the cost of integrating along a vertical line from
 * low_lim to hi_lim. The variable x represents the incomming position
 * on the bottom boundary and y represents the position on the outgoing boundary.
 * 
 * @param low_lim 
 * @param hi_lim 
 * @param cell 
 * @param x_coef 
 * @param x_const boolean determining whether x is to be treated as a variable
 * @return std::vector<ConstrainedBivariatePolynomial<2>> 
 */
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

    double sx = cell.s.x;
    double sy = cell.s.y;
    double tx = cell.t.x;
    double ty = cell.t.y;

    double alpha = angle(s1, t1);
    double beta = angle(s2, t2);
    
    if (x_const) {
        auto terms_x = solve_integral(low_lim, hi_lim,
            // polynomial
            BivariatePolynomial<1>({{{{x_coef*_cos(alpha) + s1.x-s2.x, 0}},{{0, 0}}}}),
            // a
            _cos(beta), 
            {sy, ty},
            {Polynomial<1>({sx, 0})}, {Polynomial<1>({tx, 0})}
        );

        for (auto term: terms_x)
            xterms.push_back(term);

        auto terms_y = solve_integral(low_lim, hi_lim,
            // polynomial
            BivariatePolynomial<1>({{{{x_coef * _sin(alpha) + s1.y-s2.y, 0}},{{0, 0}}}}),
            // a
            _sin(beta), 
            {sy, ty},
            {Polynomial<1>({sx, 0})}, {Polynomial<1>({tx, 0})}
        );

        for (auto term: terms_y)
            yterms.push_back(term);
    
    } else {

        auto integrad = BivariatePolynomial<1>({{{{s1.x-s2.x, 0}},{{_cos(alpha), 0}}}});
        auto a = _cos(beta);

        auto terms_x = solve_integral(low_lim, hi_lim,
            BivariatePolynomial<1>({{{{s1.x-s2.x, 0}},{{_cos(alpha), 0}}}}),
            _cos(beta), {sy, ty},
            {Polynomial<1>({sx, 0})}, {Polynomial<1>({tx, 0})}
        );

        for (auto term: terms_x)
            xterms.push_back(term);
            
        auto terms_y = solve_integral(low_lim, hi_lim,
            BivariatePolynomial<1>({{{{s1.y-s2.y, 0}},{{_sin(alpha), 0}}}}),
            _sin(beta), {sy, ty},
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

/**
 * @brief Computes a vector of constrained bivariate polynomials
 * representing the cost of integrating along the specified axis
 * from low_liw to hi_lim.
 * 
 * @param low_lim 
 * @param hi_lim 
 * @param cell 
 * @param axis 
 * @param left_constraints 
 * @param right_constraints 
 * @param y_range 
 * @return std::vector<ConstrainedBivariatePolynomial<2>> 
 */
std::vector<ConstrainedBivariatePolynomial<2>>
axis_int(BivariatePolynomial<1> low_lim, BivariatePolynomial<1> hi_lim, const Cell& cell, Line axis,
std::vector<Polynomial<1>> left_constraints, std::vector<Polynomial<1>> right_constraints, 
Interval_c y_range) {
    Point s1 = cell.s1;
    Point s2 = cell.s2;
    Point t1 = cell.t1;
    Point t2 = cell.t2;

    double sx = cell.s.x;
    double sy = cell.s.y;
    double tx = cell.t.x;
    double ty = cell.t.y;

    double alpha = angle(s1, t1);
    double beta = angle(s2, t2);

    auto gi = axis.grad_int();

    double m = gi.x;
    double b = gi.y;

    // The integral has the following form:
    // int_{l0}^{l1} | h(t, mt+b) |dt 
    // = int_{l0}^{l1} | t*_cos(alpha) - (mt+b)*_cos(beta) | + |t*_sin(alpha) - (mt+b)*_sin(alpha)|dt
    //

    // Integral corresponding to the first term in the integrand above
    auto xterms = solve_integral(
        low_lim, hi_lim,
        // polynomial
        BivariatePolynomial<1>({{{{s1.x-s2.x - b*_cos(beta), 0}},{{0, 0}}}}),
        // a
        m*_cos(beta) - _cos(alpha), 
        y_range,
        left_constraints, right_constraints
    );

    // Integral corresponding to the second term in the integrand above
    auto yterms = solve_integral(
        low_lim, hi_lim,
        // polynomial
        BivariatePolynomial<1>({{{{s1.y-s2.y - b*_sin(beta), 0}},{{0, 0}}}}),
        // a
        m*_sin(beta) - _sin(alpha), 
        y_range,
        left_constraints, right_constraints
    );

    std::vector<ConstrainedBivariatePolynomial<2>> results = std::vector<ConstrainedBivariatePolynomial<2>>();

    for (auto xterm: xterms)
        for (auto yterm: yterms)
            results.push_back(xterm + yterm);

    auto valid_results = std::vector<ConstrainedBivariatePolynomial<2>>();

    for (auto f: results)
        if (valid_constraints(f)) {
            valid_results.push_back(f.multiply(fabs(1+m)));
        }

    return valid_results;
}

std::vector<ConstrainedBivariatePolynomial<2>>
combine_steps(size_t step_count, std::vector<std::vector<ConstrainedBivariatePolynomial<2>>> steps) {
    
    auto results = std::vector<ConstrainedBivariatePolynomial<2>>();
    
    if (step_count == 2) 
        for (size_t i = 0; i < steps[0].size(); ++i)
            for (size_t j = 0; j < steps[1].size(); ++j) {
                results.push_back(steps[0][i] + steps[1][j]);

                // debugging
                // auto poly = results.back();
                // results.back().ymax = poly.y_interval.max;
                // double xmin = -1;
                // for (auto c: poly.left_constraints) {
                //     xmin = std::max(xmin, c(poly.y_interval.max));
                // }
                // results.back().xmin = xmin;
                // auto cost = poly.slice_at_y(poly.y_interval.max).polynomial(xmin);
                // results.back().test_value = cost;
            }

    else if (step_count == 3)
        for (size_t i = 0; i < steps[0].size(); ++i)
            for (size_t j = 0; j < steps[1].size(); ++j)
                for (size_t k = 0; k < steps[2].size(); ++k) {
                    auto lhs = steps[0][i] + steps[1][j];
                    results.push_back(lhs + steps[2][k]);
                    
                    // debugging
                    // auto poly = results.back();
                    // results.back().ymax = poly.y_interval.max;
                    // double xmin = -1;
                    // for (auto c: poly.left_constraints) {
                    //     xmin = std::max(xmin, c(poly.y_interval.max));
                    // }
                    // results.back().xmin = xmin;
                    // auto cost = poly.slice_at_y(poly.y_interval.max).polynomial(xmin);
                    // results.back().test_value = cost;
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

/**
 * @brief Computes a vector of constrained bivariate quadratics respresenting
 * the cost of integrating along the cell axis 
 * 
 * @param axis 
 * @param cell 
 * @return * std::vector<ConstrainedBivariatePolynomial<2>> 
 */
std::vector<ConstrainedBivariatePolynomial<2>>
bottom_right_axis_integrals(Line axis, const Cell& cell) {
    double sx = cell.s.x; //- cell.mid.x;
    double sy = cell.s.y; //- cell.mid.y;
    double tx = cell.t.x; //- cell.mid.x;
    double ty = cell.t.y; //- cell.mid.y;

    // axis given by y = m*x + b
    double m = axis.grad_int().x;
    double b = axis.grad_int().y;

    auto costs = std::vector<ConstrainedBivariatePolynomial<2>>();

    // Impose constraints for each type of path.
    // e.g., if the path moves up to the axis, the starting position
    // must be to the right/below the axis which can be expressed as
    // x >= -b/m.
    // If the path moves horizontally from the axis to the right cell boundary,
    // the outgoing point must be below the axis which can be expressed as 
    // y <= m*tx + b
    Polynomial<1> up_to_axis_xl = Polynomial<1>({-b/m, 0});
    Polynomial<1> up_to_axis_xr = Polynomial<1>({(ty-b)/m, 0});
    double right_from_axis_ymin = 0;
    double right_from_axis_ymax = m*tx + b;

    Polynomial<1> right_to_axis_xl = Polynomial<1>({0, 0});
    Polynomial<1> right_to_axis_xr = Polynomial<1>({-b/m, 0});
    double up_from_axis_ymin = m*tx + b;
    double up_from_axis_ymax = ty;

    std::vector<Polynomial<1>> right_to_axis_l_const = {};

    // up, axis, across

    auto uaa_up_terms = vertical_int(
        BivariatePolynomial<1>({{{{0, 0}},{{0, 0}}}}),
        BivariatePolynomial<1>({{{{b, 0}},{{m, 0}}}}),
        cell, 1, false
    );

    double x_int = axis.getX(sy);
    double y_int = axis.getY(tx);

    std::vector<Polynomial<1>> test = {Polynomial<1>({tx, 0}), up_to_axis_xr, Polynomial<1>({-b/m, 1/m})};

    auto uaa_axis_terms = axis_int(
        BivariatePolynomial<1>({{{{0, 0}},{{1, 0}}}}),
        BivariatePolynomial<1>({{{{-b/m, 1/m}},{{0, 0}}}}),
        cell, axis,
        {Polynomial<1>({sx, 0}), up_to_axis_xl}, 
        {Polynomial<1>({tx, 0}), up_to_axis_xr, Polynomial<1>({-b/m, 1/m})},
        {sy, sy ? std::min(ty, y_int) < sy : std::min(ty, y_int)}
    );

    auto uaa_across_terms = horizontal_int(
        BivariatePolynomial<1>({{{{-b/m, 1/m}},{{0, 0}}}}),
        BivariatePolynomial<1>({{{{tx, 0}},{{0, 0}}}}),
        cell, 1, false
    );

    auto up_axis_across = combine_steps(3, {uaa_up_terms, uaa_axis_terms, uaa_across_terms});
    for (auto& poly: up_axis_across) {
        poly.path_type = UAA;
        poly.axis = axis;
        costs.push_back(poly);
    }

    // across-axis-up

    auto aau_across_terms = horizontal_int(
        BivariatePolynomial<1>({{{{0, 0}},{{1, 0}}}}),
        BivariatePolynomial<1>({{{{(sy-b)/m, 0}},{{0, 0}}}}),
        cell, sy, true
    );

    auto aau_axis_terms = axis_int(
        BivariatePolynomial<1>({{{{(sy-b)/m, 0}},{{0, 0}}}}),
        BivariatePolynomial<1>({{{{tx, 0}},{{0, 0}}}}),
        cell, axis,
        {Polynomial<1>({sx, 0}), right_to_axis_xl}, 
        {Polynomial<1>({tx, 0}), Polynomial<1>({(ty-b)/m, 0}), right_to_axis_xr},
        {std::max(sy, m*tx+b), ty}
    );

    auto aau_up_terms = vertical_int(
        BivariatePolynomial<1>({{{{m*tx+b, 0}},{{0, 0}}}}),
        BivariatePolynomial<1>({{{{0, 1}},{{0, 0}}}}),
        cell, tx, true);

    auto across_axis_up = combine_steps(3, {aau_up_terms, aau_axis_terms, aau_across_terms});
    for (auto& poly: across_axis_up) {
        poly.path_type = AAU;
        poly.axis = axis;
        costs.push_back(poly);
    }

    // up-axis-up
    
    auto uau_up_terms_1 = vertical_int(
        BivariatePolynomial<1>({{{{0, 0}},{{0, 0}}}}),
        BivariatePolynomial<1>({{{{b, 0}},{{m, 0}}}}),
        cell, 1, false
    );

    auto uau_axis_terms = axis_int(
        BivariatePolynomial<1>({{{{0, 0}},{{1, 0}}}}),
        BivariatePolynomial<1>({{{{tx, 0}},{{0, 0}}}}),
        cell, axis,
        {Polynomial<1>({sx, 0}), Polynomial<1>({-b/m, 0}), up_to_axis_xl}, 
        {Polynomial<1>({tx, 0}), Polynomial<1>({(ty-b)/m, 0}), up_to_axis_xr},
        {std::max(sy, m*tx+b), ty}
    );

    auto uau_up_terms_2 = vertical_int(
        BivariatePolynomial<1>({{{{m*tx+b, 0}},{{0, 0}}}}),
        BivariatePolynomial<1>({{{{0, 1}},{{0, 0}}}}),
        cell, tx, true
    );

    auto up_axis_up = combine_steps(3, {uau_up_terms_1, uau_axis_terms, uau_up_terms_2});
    for (auto& poly: up_axis_up) {
        poly.path_type = UAU;
        poly.axis = axis;
        costs.push_back(poly);
    }

    // across-axis-across

    auto aaa_across_terms_1 = horizontal_int(
        BivariatePolynomial<1>({{{{0, 0}},{{1, 0}}}}),
        BivariatePolynomial<1>({{{{(sy-b)/m, 0}},{{0, 0}}}}),
        cell, sy, true
    );

    auto aaa_axis_terms = axis_int(
        BivariatePolynomial<1>({{{{(sy-b)/m, 0}},{{0, 0}}}}),
        BivariatePolynomial<1>({{{{-b/m, 1/m}},{{0, 0}}}}),
        cell, axis,
        {Polynomial<1>({sx, 0}), right_to_axis_xl}, 
        {Polynomial<1>({tx, 0}), Polynomial<1>({-b/m, 1/m}), right_to_axis_xr},
        {sy, sy ? std::min(ty, y_int) < sy : std::min(ty, y_int)}
    );

    auto aaa_across_terms_2 = horizontal_int(
        BivariatePolynomial<1>({{{{-b/m, 1/m}},{{0, 0}}}}),
        BivariatePolynomial<1>({{{{tx, 0}},{{0, 0}}}}),
        cell, 1, false
    );

    auto across_axis_across = combine_steps(3, {aaa_across_terms_1, aaa_axis_terms, aaa_across_terms_2});
    for (auto& poly: across_axis_across) {
        poly.path_type = AAA;
        poly.axis = axis;
        costs.push_back(poly);
    }


    return costs;
}

/**
 * @brief Computes a vector of constrained bivariate quadratics
 * representing the costs of integrating along paths from the bottom
 * to the right boundary of a cell.
 * 
 * @param cell 
 * @return std::vector<ConstrainedBivariatePolynomial<2>> 
 */
std::vector<ConstrainedBivariatePolynomial<2>>
bottom_to_right_costs_2D(const Cell& cell) {
    double sx = cell.s.x;
    double sy = cell.s.y; 
    double tx = cell.t.x;
    double ty = cell.t.y;

    std::vector<ConstrainedBivariatePolynomial<2>> costs = std::vector<ConstrainedBivariatePolynomial<2>>();

    auto axes = get_axis(cell);

    double alpha = angle(cell.s1, cell.t1);
    double beta = angle(cell.s2, cell.t2);

    Point s1 = cell.s1;
    Point t1 = cell.t1;
    Point s2 = cell.s2;
    Point t2 = cell.t2;

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

    for (auto poly: up_across) {
        poly.path_type = UA;
        costs.push_back(poly);
    }

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


    for (auto poly: across_up) {
        poly.path_type = AU;
        costs.push_back(poly);
    }

    if (axes.size() == 1) {
        Line axis = axes[0];
        auto axis_costs = bottom_right_axis_integrals(axis, cell);
        for (auto poly: axis_costs)
            costs.push_back(poly);
    }

    for (auto& cost: costs)
        cost.boundaries = BR;


    return costs;
}

/**
 * @brief Computes a vector of constrained bivariate quadratics representing
 * the cost of integrating along a vertical segment where the variables x and y
 * represent positions on the segment corresponding to the bottom/top cell boundaries.
 * 
 * @param low_lim 
 * @param hi_lim 
 * @param cell 
 * @param fixed_var 
 * @return std::vector<ConstrainedBivariatePolynomial<2>> 
 */
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
            BivariatePolynomial<1>({{{{s1.x-s2.x, _cos(alpha)}},{{0, 0}}}}),
            _cos(beta), {sy, ty},
            {Polynomial<1>({sx, 0})}, {Polynomial<1>({tx, 0})} 
        );
        y_terms = solve_integral(low_lim, hi_lim,
            BivariatePolynomial<1>({{{{s1.y-s2.y, _sin(alpha)}},{{0, 0}}}}),
            _sin(beta), {sy, ty},
            {Polynomial<1>({sx, 0})}, {Polynomial<1>({tx, 0})} 
        );
    } else if (fixed_var == "x") {
        x_terms = solve_integral(low_lim, hi_lim,
            BivariatePolynomial<1>({{{{s1.x-s2.x, 0}},{{_cos(alpha), 0}}}}),
            _cos(beta), {sy, ty},
            {Polynomial<1>({sx, 0})}, {Polynomial<1>({tx, 0})} 
        );
        y_terms = solve_integral(low_lim, hi_lim,
            BivariatePolynomial<1>({{{{s1.y-s2.y, 0}},{{_sin(alpha), 0}}}}),
            _sin(beta), {sy, ty},
            {Polynomial<1>({sx, 0})}, {Polynomial<1>({tx, 0})} 
        );
    }


    for (auto xterm: x_terms)
        for (auto yterm: y_terms)
            results.push_back(xterm + yterm);

    auto valid_results = std::vector<ConstrainedBivariatePolynomial<2>>();

    for (auto f: results)
        if (valid_constraints(f)) {
            valid_results.push_back(f);
        }

    return valid_results;
}

/**
 * @brief Computes constrained bivariate quadratics
 * representing the cost of integrating along a horizontal
 * line where both variables x and y represent positions on the
 * the segment corresponding to the bottom/top cell boundaries.
 * 
 * @param low_lim 
 * @param hi_lim 
 * @param cell 
 * @param height 
 * @return std::vector<ConstrainedBivariatePolynomial<2>> 
 */
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
        // polynomial
        BivariatePolynomial<1>({{{{s2.x-s1.x + height*_cos(beta), 0}},{{0, 0}}}}),
        // a
        _cos(alpha), 
        {sy, ty},
        {Polynomial<1>({sx, 0})}, {Polynomial<1>({tx, 0})}
    );

    auto y_terms = solve_integral(
        low_lim, hi_lim,
        // polynomial
        BivariatePolynomial<1>({{{{s2.y-s1.y + height*_sin(beta), 0}},{{0, 0}}}}),
        // a
        _sin(alpha), 
        {sy, ty},
        {Polynomial<1>({sx, 0})}, {Polynomial<1>({tx, 0})}
    );

    for (auto xterm: x_terms)
        for (auto yterm: y_terms)
            results.push_back(xterm + yterm);

    auto valid_results = std::vector<ConstrainedBivariatePolynomial<2>>();

    for (auto f: results)
        if (valid_constraints(f)) {
            valid_results.push_back(f);
        }

    return valid_results;
}

/**
 * @brief Computes a vector of constrained bivariate quadratics
 * representing the cost of integrating along the cell axis where
 * both x and y variables are positions on the segment corresponding to the bottom/top
 * cell boundaries
 * 
 * @param axis 
 * @param cell 
 * @return std::vector<ConstrainedBivariatePolynomial<2>> 
 */
std::vector<ConstrainedBivariatePolynomial<2>>
bottom_top_axis_integrals(Line axis, const Cell& cell) {

    auto sx = 0.;
    auto tx = cell.len1;
    auto sy = 0.;
    auto ty = cell.len2;

    double alpha = angle(cell.s1, cell.t1);
    double beta = angle(cell.s2, cell.t2);

    Point s1 = cell.s1;
    Point t1 = cell.t1;
    Point s2 = cell.s2;
    Point t2 = cell.t2;

    // The axis is given by y = mx + b
    double m = axis.grad_int().x;
    double b = axis.grad_int().y;

    auto results = std::vector<ConstrainedBivariatePolynomial<2>>();

    auto bottom_int = axis.getX(0);
    auto top_int = axis.getX(cell.len2);

    // Constraints on x which need to be satisfied for each path type
    // eg: up_to_axis_xl = left constraint on x for a path which moves up to the axis
    // from the bottom of the cell. 
    // For such a path, the starting point, (x, 0), must
    // be to the right of the axis which can be expressed as -b/m <= x.
    Polynomial<1> up_to_axis_xl = Polynomial<1>({-b/m, 0});
    Polynomial<1> up_to_axis_xr = Polynomial<1>({(ty-b)/m, 0});
    Polynomial<1> right_to_axis_xl = Polynomial<1>({0, 0});
    Polynomial<1> right_to_axis_xr = Polynomial<1>({-b/m, 0});
    
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
        {Polynomial<1>({sx, 0}), Polynomial<1>({bottom_int}), up_to_axis_xl}, 
        {Polynomial<1>({tx, 0}), up_to_axis_xr},
        {sx, std::min(tx, top_int)}
    );

    auto uau_u2_terms = vertical_int_b2t(
        BivariatePolynomial<1>({{{{b, m}},{{0, 0}}}}),
        BivariatePolynomial<1>({{{{cell.len2 , 0}},{{0, 0}}}}),
        cell, "y"
    );

    auto up_axis_up = combine_steps(3, {uau_u1_terms, uau_axis_terms, uau_u2_terms});
    for (auto poly: up_axis_up) {
        poly.path_type = UAU;
        poly.axis = axis;
        results.push_back(poly);
    }

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
        {Polynomial<1>({sx, 0}), right_to_axis_xl}, 
        {Polynomial<1>({tx, 0}), Polynomial<1>({bottom_int, 0}), right_to_axis_xr},
        {sx, std::min(tx, top_int)}
    );

    auto aau_up_terms = vertical_int_b2t(
        BivariatePolynomial<1>({{{{b, m}},{{0, 0}}}}),
        BivariatePolynomial<1>({{{{cell.len2 , 0}},{{0, 0}}}}),
        cell, "y"
    );

    auto across_axis_up = combine_steps(3, {aau_across_terms, aau_axis_terms, aau_up_terms});
    for (auto poly: across_axis_up) {
        poly.path_type = AAU;
        poly.axis = axis;
        results.push_back(poly);
    }

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
        {Polynomial<1>({sx, 0}), Polynomial<1>({bottom_int, 0}), up_to_axis_xl}, 
        {Polynomial<1>({tx, 0}), up_to_axis_xr},
        {std::max(sx, top_int), tx}
    );

    auto uaa_across_terms = horizontal_int_b2t(
        BivariatePolynomial<1>({{{{top_int, 0}},{{0, 0}}}}),
        BivariatePolynomial<1>({{{{0, 1}},{{0, 0}}}}),
        cell, cell.len2
    );

    auto up_axis_across = combine_steps(3, {uaa_up_terms, uaa_axis_terms, uaa_across_terms});
    for (auto poly: up_axis_across) {
        poly.path_type = UAA;
        poly.axis = axis;
        results.push_back(poly);
    }

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
        {Polynomial<1>({sx, 0}), right_to_axis_xl}, 
        {Polynomial<1>({tx, 0}), Polynomial<1>({bottom_int, 0}), right_to_axis_xr},
        {std::max(sx, top_int), tx}
    );

    auto aaa_a2_terms = horizontal_int_b2t(
        BivariatePolynomial<1>({{{{top_int, 0}},{{0, 0}}}}),
        BivariatePolynomial<1>({{{{0, 1}},{{0, 0}}}}),
        cell, cell.len2
    );

    auto across_axis_across = combine_steps(3, {aaa_a1_terms, aaa_axis_terms, aaa_a2_terms});
    for (auto poly: across_axis_across) {
        poly.path_type = AAA;
        poly.axis = axis;
        results.push_back(poly);
    }

    auto valid_results = std::vector<ConstrainedBivariatePolynomial<2>>();

    for (auto f: results)
        if (valid_constraints(f))
            valid_results.push_back(f);

    for (auto f: valid_results)
        f.boundaries = BT;

    return valid_results;
}

std::vector<ConstrainedBivariatePolynomial<2>>
bottom_to_top_costs_2D(const Cell& cell) {
    double sx = cell.s.x;
    double sy = cell.s.x;
    double tx = cell.t.x;
    double ty = cell.t.x;

    std::vector<ConstrainedBivariatePolynomial<2>> costs = std::vector<ConstrainedBivariatePolynomial<2>>();

    auto axes = get_axis(cell);

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

    auto up_across = combine_steps(2, {ua_up_terms, ua_across_terms});
    for (auto poly: up_across) {
        poly.path_type = UA;
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
        poly.path_type = AU;
        costs.push_back(poly);
    }

    if (axes.size() == 1) {
        Line axis = axes[0];
        auto axis_costs = bottom_top_axis_integrals(axis, cell);
        for (auto poly: axis_costs) {
            costs.push_back(poly);
        }
    }
    
    for (auto& cost: costs)
        cost.boundaries = BT;

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
                    total_cost.right_constraints,
                    cell_cost
                );

                min_pieces.insert(min_pieces.end(), min_total_cost.pieces.begin(), min_total_cost.pieces.end());
            }
            ++pieces_it;
        }


        const auto result = FAST_LOWER_ENVELOPE ? fast_lower_envelope_v2(min_pieces) 
        : naive_lower_envelope(min_pieces);


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
    double sx = cell.s.x; 
    double sy = cell.s.y; 
    double tx = cell.t.x;
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
    
    const auto final_result = FAST_LOWER_ENVELOPE ? fast_lower_envelope_v2(results.pieces) 
        : naive_lower_envelope(results.pieces);

    return final_result;
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
    auto costs = bottom_to_right_costs_2D(cell);
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
    return bottom_to_top_costs_2D(cell);
}

#endif //TRAJECTORY_CLUSTERING_2D_L1_L1_H