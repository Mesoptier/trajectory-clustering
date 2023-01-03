#include "integral_solver.h"

std::vector<Polynomial<1>> x_constraints(
    BivariatePolynomial<1> lhs, BivariatePolynomial<1> rhs, 
    std::string ineq, std::string side) {
        auto lx = lhs.coefficients[1][0];
        auto ly = lhs.coefficients[0][1];
        auto lc = lhs.coefficients[0][0];
        auto rx = rhs.coefficients[1][0];
        auto ry = rhs.coefficients[0][1];
        auto rc = rhs.coefficients[0][0];

        if (side == "left") {
            if ((ineq == ">" && lx < rx) || 
                (ineq == "<" && lx > rx))
                return {};
        }   else if (side == "right") {
            if ((ineq == ">" && lx > rx) || 
                (ineq == "<" && lx < rx))
                return {};
        }
        
        if (approx_equal(lx, rx))
            return {};

        return {Polynomial<1>({{
            (rc - lc) / (lx - rx), (ry - ly) / (lx - rx)
        }})};
}


std::vector<double> y_constraints(
    BivariatePolynomial<1> lhs, BivariatePolynomial<1> rhs, 
    std::string ineq, std::string side) {
        auto lx = lhs.coefficients[1][0];
        auto ly = lhs.coefficients[0][1];
        auto lc = lhs.coefficients[0][0];
        auto rx = rhs.coefficients[1][0];
        auto ry = rhs.coefficients[0][1];
        auto rc = rhs.coefficients[0][0];

        if (!approx_equal(lx, rx))
            return {};

        if (side == "upper") {
            if ((ineq == "<" && ly < ry) || 
                (ineq == ">" && ly > ry))
                    return {};
            if (approx_equal(ly, ry)) {
                if (ineq == "<" && (lc - rc < 0 || approx_equal(lc, rc)))
                    return {};
                else if (ineq == ">" && (lc - rc > 0 || approx_equal(lc, rc)))
                    return {};
                else
                    return {-1};
            }
        } else if (side == "lower") {
            if ((ineq == "<" && ly > ry) || 
                (ineq == ">" && ly < ry) ||
                approx_equal(ly, ry))
                    return {};
        }


        return {(rc - lc) / (ly - ry)};    
}

void update_constraints(BivariatePolynomial<1> lhs, BivariatePolynomial<1> rhs, std::string ineq,
    std::vector<Polynomial<1>>& left_constraints, std::vector<Polynomial<1>>& right_constraints,
    double& y_lower_bound, double& y_upper_bound) {
    
    bool constraint_found = false;

    if (!constraint_found) {
        for (auto left: x_constraints(lhs, rhs, ineq, "left")) {
            left_constraints.push_back(left);
            constraint_found = true;
        }
    } 
    if (!constraint_found) {
        for (auto right: x_constraints(lhs, rhs, ineq, "right")) {
            right_constraints.push_back(right);
            constraint_found = true;
        }
    } 
    if (!constraint_found) {
        for (auto upper: y_constraints(lhs, rhs, ineq, "upper")) {
            y_upper_bound = std::min(y_upper_bound, upper);
            constraint_found = true;
        }
    }
    if (!constraint_found) {
        for (auto lower: y_constraints(lhs, rhs, ineq, "lower")) {
            y_lower_bound = std::max(y_lower_bound, lower);
            constraint_found = true;
        }
    }
}

/**
 * Explicity solves an integral and returns a vector
 * of constrained bivariate quadratic functions
*/
std::vector<ConstrainedBivariatePolynomial<2>>
solve_integral(BivariatePolynomial<1> low_lim, 
BivariatePolynomial<1> hi_lim,
BivariatePolynomial<1> integrand, double t_coeff,
Interval_c y_range, std::vector<Polynomial<1>> left_constraints, std::vector<Polynomial<1>> right_constraints) {
    auto output = std::vector<ConstrainedBivariatePolynomial<2>>();
        /*
        This function solves an integral of the form
        int_{low_lim}^{hi_lim}|integrad - a*t|dt
        where low_lim, hi_lim and integrad are linear
        functions in two variables.
        */
    


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

        /* 
        We use the following notation to compute the integral:
        low_lim = P_1(x, y) = p1x*x + p1y*y + p1c
        hi_lim = P_2(x, y) = p2x*x + p2y*y + p2c
        integrad = P_3(x, y) = p3x*x + p3y*y + p3c
        The integral then becomes int_{P_1}^{P_2}|P_3 - a*t|dt
        */
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

        BivariatePolynomial<1> P1 = BivariatePolynomial<1>({{
            {p1c, p1y}, {p1x, 0}
        }});
        BivariatePolynomial<1> P2 = BivariatePolynomial<1>({{
            {p2c, p2y}, {p2x, 0}
        }});
        BivariatePolynomial<1> P3 = BivariatePolynomial<1>({{
            {p3c, p3y}, {p3x, 0}
        }});

        // Case 1: P_3(x,y) - a*P_1(x,y) >= 0 and P_3(x,y) - a*P_2(x,y) >= 0
        // In this case, we compute the integral by simply dropping the absolute value sign.
        auto case_1_poly = BivariatePolynomial<2>({{
            {
                // x^0 terms
                {
                    // constant 
                    a*pow(p1c, 2)/2 - a*pow(p2c, 2)/2 - p1c*p3c + p2c*p3c,
                    // y 
                    a*p1c*p1y - a*p2c*p2y - p1c*p3y - p1y*p3c + p2c*p3y + p2y*p3c,
                    // y^2 
                    a*pow(p1y, 2)/2 - a*pow(p2y, 2)/2 - p1y*p3y + p2y*p3y
                }
            },
            {
                // x^1 terms
                {
                    //x
                    a*p1c*p1x - a*p2c*p2x - p1c*p3x - p1x*p3c + p2c*p3x + p2x*p3c,
                    // xy
                    a*p1x*p1y - a*p2x*p2y - p1x*p3y - p1y*p3x + p2x*p3y + p2y*p3x,
                    //xy^2
                    0
                }
            },
            {
                // x^2 terms
                {
                    // x^2
                    a*pow(p1x, 2)/2 - a*pow(p2x, 2)/2 - p1x*p3x + p2x*p3x,
                    // x^2y
                    0,
                    // x^2y^2
                    0
                }
            }
        }});
        
        auto case_1_right_constraints = std::vector<Polynomial<1>>();
        auto case_1_left_constraints = std::vector<Polynomial<1>>();

        double y_upper_bound = y_range.max;
        double y_lower_bound = y_range.min;

        // Impose constraints on x and y to ensure P_3(x,y) - a*P_1(x,y) >= 0 and P_3(x,y) - a*P_2(x,y) >= 0
        // is satisfied.
        // ineq1: p3x*x + p3y*y + p3c -a*p1x*x - a*p1y*y - a*p1c >= 0
        //        (p3x - a*p1x)*x >= (a*p1 - p3y)*y + (a*p1c - p3c)
        // ineq2: p3x*x + p3y*y + p3c -a*p2x*x - a*p2y*y - a*p2c >= 0
        //        (p3x - a*p2x)*x >= (a*p2 - p3y)*y + (a*p2c - p3c)

        // if (approx_equal(p3x, a*p2x)) {
        //     if (a*p2y > p3y && !approx_equal(a*p2y, p3y)) {
        //         y_upper_bound = std::min(y_upper_bound, (p3c-a*p2c) / (a*p2y-p3y));
        //     }
        //     else if (a*p2y < p3y && !approx_equal(a*p2y, p3y))
        //         y_lower_bound = std::max(y_lower_bound, (p3c-a*p2c) / (a*p2y-p3y));
        //     else if (approx_equal(a*p2y, p3y)) {
        //         if (!(p3c-a*p2c > 0) && 
        //         !approx_equal(p3c, a*p2c)) {
        //             y_upper_bound = -1;
        //         }
        //     }
        // } 
        // else if (p3x > a*p2x) {
        //     case_1_left_constraints.push_back(Polynomial<1>({
        //         (a*p2c - p3c) / (p3x-a*p2x), 
        //         (a*p2y-p3y) / (p3x-a*p2x) 
        //     }));
        // } else if (p3x < a*p2x) {
        //     case_1_right_constraints.push_back(Polynomial<1>({
        //         (a*p2c-p3c) / (p3x-a*p2x), 
        //         (a*p2y-p3y) / (p3x-a*p2x) 
        //     }));
        // } 

        update_constraints(P3, P2*a, ">", case_1_left_constraints, case_1_right_constraints,
        y_lower_bound, y_upper_bound);
        

        // if (approx_equal(p3x, a*p1x)) {
        //     if (a*p1y > p3y) {
        //         y_upper_bound = std::min(y_upper_bound, (p3c-a*p1c) / (a*p1y-p3y));
        //     }
        //     else if (a*p1y < p3y)
        //         y_lower_bound = std::max(y_lower_bound, (p3c-a*p1c) / (a*p1y-p3y));
        //     else if (approx_equal(a*p1y, p3y)) {
        //         if (!(p3c-a*p1c > 0 || approx_equal(p3c, a*p1c))) {
        //             y_upper_bound = -1;
        //         }
        //     }
        // } else if (p3x > a*p1x) {
        //     case_1_left_constraints.push_back(Polynomial<1>({
        //         (a*p1c-p3c) / (p3x-a*p1x), 
        //         (a*p1y-p3y) / (p3x-a*p1x) 
        //     }));
        // } else if (p3x < a*p1x) {
        //     case_1_right_constraints.push_back(Polynomial<1>({
        //         (a*p1c-p3c) / (p3x-a*p1x), 
        //         (a*p1y-p3y) / (p3x-a*p1x) 
        //     }));
        // }

        update_constraints(P3, P1*a, ">", case_1_left_constraints, case_1_right_constraints,
        y_lower_bound, y_upper_bound);

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


        output.back().history = "1";

        // Case 2: P_3(x, y) -a*P_1(x, y) <= 0 && P_3(x, y) - a*P_2(x, y) <= 0
        // In this case, we compute the integral by simply dropping the absolute value sign
        // and multiplying by -1.
        auto case_2_poly = BivariatePolynomial<2>({{
            {   
                // x^0 terms
                {
                    // consant
                    -a*pow(p1c, 2)/2 + a*pow(p2c, 2)/2 + p1c*p3c - p2c*p3c,
                    // y
                    -a*p1c*p1y + a*p2c*p2y + p1c*p3y + p1y*p3c - p2c*p3y - p2y*p3c,
                    // y^2
                    -a*pow(p1y, 2)/2 + a*pow(p2y, 2)/2 + p1y*p3y - p2y*p3y
                }
            },
            {   
                // x^1 terms
                {
                    // x
                    -a*p1c*p1x + a*p2c*p2x + p1c*p3x + p1x*p3c - p2c*p3x - p2x*p3c,
                    // xy
                    -a*p1x*p1y + a*p2x*p2y + p1x*p3y + p1y*p3x - p2x*p3y - p2y*p3x,
                    // xy^2
                    0
                }
            },
            {   
                // x^2 terms
                {
                    // x^2
                    -a*pow(p1x, 2)/2 + a*pow(p2x, 2)/2 + p1x*p3x - p2x*p3x,
                    // x^2y
                    0,
                    // x^2y^2
                    0
                }
            }        
        }});
        
        auto case_2_right_constraints = std::vector<Polynomial<1>>();
        auto case_2_left_constraints = std::vector<Polynomial<1>>();

        y_upper_bound = y_range.max;
        y_lower_bound = y_range.min;

        // Impose constraints to ensure P_3(x, y) -a*P_1(x, y) <= 0 && P_3(x, y) - a*P_2(x, y) <= 0
        // if (approx_equal(p3x, a*p2x)) {
        //     if (a*p2y > p3y)
        //         y_lower_bound = std::max(y_lower_bound, (a*p2c-p3c) / (p3y-a*p2y));
        //     else if (a*p2y < p3y)
        //         y_upper_bound = std::min(y_upper_bound, (a*p2c-p3c) / (p3y-a*p2y));
        //     else if (approx_equal(a*p2y, p3y)) {
        //         if (!(p3c-a*p2c < 0) && !(approx_equal(p3c, a*p2c))) {
        //             y_upper_bound = -1;
        //         }
        //     }
        // } else if (p3x > a*p2x) {
        //     case_2_right_constraints.push_back(Polynomial<1>({
        //         (a*p2c-p3c) / (p3x-a*p2x), 
        //         (a*p2y-p3y) / (p3x-a*p2x) 
        //     }));
        // } else if (p3x < a*p2x) {
        //     case_2_left_constraints.push_back(Polynomial<1>({
        //         (a*p2c-p3c) / (p3x-a*p2x), 
        //         (a*p2y-p3y) / (p3x-a*p2x) 
        //     }));
        // }

        update_constraints(P3, P2*a, "<", case_2_left_constraints, case_2_right_constraints,
        y_lower_bound, y_upper_bound);

        // if (approx_equal(p3x, a*p1x)) {
        //     if (a*p1y > p3y)
        //         y_lower_bound = std::max(y_lower_bound, (p3c-a*p1c) / (a*p1y-p3y));
        //     else if (a*p1y < p3y)
        //         y_upper_bound = std::min(y_upper_bound, (p3c-a*p1c) / (a*p1y-p3y));
        //     else if (approx_equal(a*p1y, p3y)) {
        //         if (!(p3c-a*p1c < 0) && !(approx_equal(p3c, a*p1c))) {
        //             y_upper_bound = -1;
        //         }
        //     }
        // } else if (p3x > a*p1x) {
        //     case_2_right_constraints.push_back(Polynomial<1>({
        //         (a*p1c-p3c) / (p3x-a*p1x), 
        //         (a*p1y-p3y) / (p3x-a*p1x) 
        //     }));
        // } else if (p3x < a*p1x) {
        //     case_2_left_constraints.push_back(Polynomial<1>({
        //         (a*p1c-p3c) / (p3x-a*p1x), 
        //         (a*p1y-p3y) / (p3x-a*p1x) 
        //     }));
        // }

        update_constraints(P3, P1*a, "<", case_2_left_constraints, case_2_right_constraints,
        y_lower_bound, y_upper_bound);

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
        
        output.back().history = "2";

        //Case 3a: P_3(x, y) - a*P_3(x, y) >= 0 and P_3(x, y) - a*P_1(x, y) <= 0
        auto case_3a_poly = BivariatePolynomial<2>({{
            {
                // x^0 terms
                {
                    // constant
                    (pow(a, 2)*pow(p1c, 2) + pow(a, 2)*pow(p2c, 2) - 2*a*p1c*p3c - 2*a*p2c*p3c + 2*pow(p3c, 2))/(2*a),
                    // y
                    (pow(a, 2)*p1c*p1y + pow(a, 2)*p2c*p2y - a*p1c*p3y - a*p1y*p3c - a*p2c*p3y - a*p2y*p3c + 2*p3c*p3y)/a,
                    // y^2
                    (pow(a, 2)*pow(p1y, 2) + pow(a, 2)*pow(p2y, 2) - 2*a*p1y*p3y - 2*a*p2y*p3y + 2*pow(p3y, 2))/(2*a)
                }
            },
            {
                // x^1 terms
                {
                    // x
                    (pow(a, 2)*p1c*p1x + pow(a, 2)*p2c*p2x - a*p1c*p3x - a*p1x*p3c - a*p2c*p3x - a*p2x*p3c + 2*p3c*p3x)/a,
                    // xy
                    (pow(a, 2)*p1x*p1y + pow(a, 2)*p2x*p2y - a*p1x*p3y - a*p1y*p3x - a*p2x*p3y - a*p2y*p3x + 2*p3x*p3y)/a,
                    // xy^2
                    0
                }
            },
            {
                // x^2 terms
                {
                    // x^2
                    (pow(a, 2)*pow(p1x, 2) + pow(a, 2)*pow(p2x, 2) - 2*a*p1x*p3x - 2*a*p2x*p3x + 2*pow(p3x, 2))/(2*a),
                    // x^2y
                    0,
                    // x^2y^2
                    0
                }
            }
        }});
        
        auto case_3a_right_constraints = std::vector<Polynomial<1>>();
        auto case_3a_left_constraints = std::vector<Polynomial<1>>();

        y_upper_bound = y_range.max;
        y_lower_bound = y_range.min;

        // if (approx_equal(p3x, p2x)) {
        //     if (a*p2y > p3y)
        //         y_lower_bound = std::max(y_lower_bound, (p3c-a*p2c) / (a*p2y-p3y));
        //     else if (a*p2y < p3y)
        //         y_upper_bound = std::min(y_upper_bound, (p3c-a*p2c) / (a*p2y-p3y));
        //     else if (approx_equal(a*p2y, p3y)) {
        //         if (!(p3c-a*p2c < 0) && !(approx_equal(p3c, a*p2c))) {
        //             y_upper_bound = -1;
        //         }
        //     }
        // } else if (p3x > a*p2x) {
        //     case_3a_right_constraints.push_back(Polynomial<1>({
        //         (a*p2c-p3c) / (p3x-a*p2x), 
        //         (a*p2y-p3y) / (p3x-a*p2x) 
        //     }));
        // } else if (p3x < a*p2x) {
        //     case_3a_left_constraints.push_back(Polynomial<1>({
        //         (a*p2c-p3c) / (p3x-a*p2x), 
        //         (a*p2y-p3y) / (p3x-a*p2x) 
        //     }));
        // }

        update_constraints(P3, P2*a, "<", case_3a_left_constraints, case_3a_right_constraints,
        y_lower_bound, y_upper_bound);

        // if (approx_equal(p3x, a*p1x)) {
        //     if (a*p1y > p3y)
        //         y_upper_bound = std::min(y_upper_bound, (p3c-a*p1c) / (a*p1y-p3y));
        //     else if (a*p1y < p3y)
        //         y_lower_bound = std::max(y_lower_bound, (p3c-a*p1c) / (a*p1y-p3y));
        //     else if (approx_equal(a*p1y, p3y)) {
        //         if (!(p3c-a*p1c > 0) && !(approx_equal(p3c, a*p1c))) {
        //             y_upper_bound = -1;
        //         }
        //     }
        // } else if (p3x > a*p1x) {
        //     case_3a_left_constraints.push_back(Polynomial<1>({
        //         (a*p1c-p3c) / (p3x-a*p1x), 
        //         (a*p1y-p3y) / (p3x-a*p1x) 
        //     }));
        // } else if (p3x < a*p1x) {
        //     case_3a_right_constraints.push_back(Polynomial<1>({
        //         (a*p1c-p3c) / (p3x-a*p1x), 
        //         (a*p1y-p3y) / (p3x-a*p1x) 
        //     }));
        // }

        update_constraints(P3, P1*a, ">", case_3a_left_constraints, case_3a_right_constraints,
        y_lower_bound, y_upper_bound);

        for (auto constraint: left_constraints) 
            case_3a_left_constraints.push_back(constraint);

         for (auto constraint: right_constraints) 
            case_3a_right_constraints.push_back(constraint);

        if (!approx_zero(a)) {
            output.push_back(ConstrainedBivariatePolynomial<2>{
                case_3a_poly,
                {y_lower_bound, y_upper_bound},
                case_3a_left_constraints,
                case_3a_right_constraints
            });

            output.back().history = "3a";
        }


        //case 3b
        auto case_3b_poly = BivariatePolynomial<2>({{
            {
                // x^0 terms
                {
                    // constant
                    (-pow(a, 2)*pow(p1c, 2) - pow(a, 2)*pow(p2c, 2) + 2*a*p1c*p3c + 2*a*p2c*p3c - 2*pow(p3c, 2))/(2*a),
                    // y
                    (-pow(a, 2)*p1c*p1y - pow(a, 2)*p2c*p2y + a*p1c*p3y + a*p1y*p3c + a*p2c*p3y + a*p2y*p3c - 2*p3c*p3y)/a,
                    // y^2
                    (-pow(a, 2)*pow(p1y, 2) - pow(a, 2)*pow(p2y, 2) + 2*a*p1y*p3y + 2*a*p2y*p3y - 2*pow(p3y, 2))/(2*a)
                }
            },
            {
                // x^1 terms
                {
                    // x
                    (-pow(a, 2)*p1c*p1x - pow(a, 2)*p2c*p2x + a*p1c*p3x + a*p1x*p3c + a*p2c*p3x + a*p2x*p3c - 2*p3c*p3x)/a,
                    // xy
                    (-pow(a, 2)*p1x*p1y - pow(a, 2)*p2x*p2y + a*p1x*p3y + a*p1y*p3x + a*p2x*p3y + a*p2y*p3x - 2*p3x*p3y)/a,
                    // xy^2
                    0
                }
            },
            {
                // x^2 terms
                {
                    // x^2
                    (-pow(a, 2)*pow(p1x, 2) - pow(a, 2)*pow(p2x, 2) + 2*a*p1x*p3x + 2*a*p2x*p3x - 2*pow(p3x, 2))/(2*a),
                    // x^2y
                    0,
                    // x^2y^2
                    0
                }
            }
        }});

        auto case_3b_right_constraints = std::vector<Polynomial<1>>();
        auto case_3b_left_constraints = std::vector<Polynomial<1>>();

        y_upper_bound = y_range.max;
        y_lower_bound = y_range.min;

        // if (approx_equal(p3x, a*p2x)) {
        //     if (a*p2y > p3y)
        //         y_upper_bound = std::min(y_upper_bound, (p3c-a*p2c) / (a*p2y-p3y));
        //     else if (a*p2y < p3y)
        //         y_lower_bound = std::max(y_lower_bound, (p3c-a*p2c) / (a*p2y-p3y));
        //     else if (approx_equal(a*p2y, p3y)) {
        //         if (!(p3c-a*p2c > 0) && !(approx_equal(p3c, a*p2c))) {
        //             y_upper_bound = -1;
        //         }
        //     }    
        // } else if (p3x > a*p2x) {
        //     case_3b_left_constraints.push_back(Polynomial<1>({
        //         (a*p2c-p3c) / (p3x-a*p2x), 
        //         (a*p2y-p3y) / (p3x-a*p2x) 
        //     }));
        // } else if (p3x < p2x) {
        //     case_3b_right_constraints.push_back(Polynomial<1>({
        //         (a*p2c-p3c) / (p3x-a*p2x), 
        //         (a*p2y-p3y) / (p3x-a*p2x) 
        //     }));
        // }

        update_constraints(P3, P2*a, ">", case_3b_left_constraints, case_3b_right_constraints,
        y_lower_bound, y_upper_bound);

        // if (approx_equal(p3x, a*p1x)) {
        //     if (a*p1y > p3y)
        //         y_lower_bound = std::max(y_lower_bound, (p3c-a*p1c) / (a*p1y-p3y));
        //     else if (a*p1y < p3y)
        //         y_upper_bound = std::min(y_upper_bound, (p3c-a*p1c) / (a*p1y-p3y));
        //     else if (approx_equal(a*p1y, p3y)) {
        //         if (!(p3c-a*p1c < 0) && !(approx_equal(p3c, a*p1c))) {
        //             y_upper_bound = -1;
        //         }
        //     }
        // } else if (p3x > a*p1x) {
        //     case_3b_right_constraints.push_back(Polynomial<1>({
        //         (a*p1c-p3c) / (p3x-a*p1x), 
        //         (a*p1y-p3y) / (p3x-a*p1x) 
        //     }));
        // } else if (p3x < a*p1x) {
        //     case_3b_left_constraints.push_back(Polynomial<1>({
        //         (a*p1c-p3c) / (p3x-a*p1x), 
        //         (a*p1y-p3y) / (p3x-a*p1x) 
        //     }));
        // }

        update_constraints(P3, P1*a, "<", case_3b_left_constraints, case_3b_right_constraints,
        y_lower_bound, y_upper_bound);

        for (auto constraint: left_constraints) 
            case_3b_left_constraints.push_back(constraint);

         for (auto constraint: right_constraints) 
            case_3b_right_constraints.push_back(constraint);

        if (!approx_zero(a)) {
            output.push_back(ConstrainedBivariatePolynomial<2>{
                case_3b_poly,
                {y_lower_bound, y_upper_bound},
                case_3b_left_constraints,
                case_3b_right_constraints
            });

            output.back().history = "3b";

        }

    auto final_output = std::vector<ConstrainedBivariatePolynomial<2>>();
        

    for (auto& poly: output)
        if (utils::valid_constraints(poly)) {
            poly.clean_constraints();
            final_output.push_back(poly);
        }

    return final_output;
}
