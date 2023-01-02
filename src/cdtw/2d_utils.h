namespace utils {
/**
 * Get the angle (measured counter-clockwise) made by the segment with 
 * the horizontal line through the source
*/
long double angle(Point s, Point t) {
    double PI = 3.14159265359;

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

/**
 * @brief return the intersection point of a and b.
 * 
 * @param a 
 * @param b 
 * @return Point 
 */
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

/**
 * @brief Check that a, b and c define a non-empty region.
 * 
 * @param a 
 * @param b 
 * @param c 
 * @return true 
 * @return false 
 */
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
            return false;
        }
    }

    return true;
}

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
        if (valid_constraints(poly)) {
            poly.clean_constraints();
            final_output.push_back(poly);
        }

    return final_output;
}


bool contains_midpoint(const Cell& cell) {
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

/**
 * Determines if axis is monote
 * @param axis
 * @return true if axis is monotone, false otherwise
 * */
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

/**
 * @brief 
 * 
 * @param cell 
 * @param p 
 * @return distance_t 
 */
distance_t l1_h_function(const Cell& cell, Point p) {
    Point vec1 = (cell.t1 - cell.s1) / (cell.t1.dist(cell.s1));
    Point vec2 = (cell.t2 - cell.s2) / (cell.t2.dist(cell.s2));
    Point p1 = cell.s1 + vec1 * p.x;
    Point p2 = cell.s2 + vec2 * p.y;
    return abs(p1.x - p2.x) + abs(p1.y - p2.y);
}

/*
* Returns the proper axis of the arugument cell.
* v1 and v2 must be the axes of the cell.
*/
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

/**
* Get axes for given cell
* @param cell
* @return vector of lines corresponding to the valleys which
* intersect the cell
*/
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