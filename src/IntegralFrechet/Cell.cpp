#include "IntegralFrechet/Cell.h"
#include "IntegralFrechet/metrics/include.h"

std::pair<Point, Point> Cell::interpolate_at(Point const& p) const {
    return {
        ImplicitEdge::interpolate_at(s1, t1, p.x),
        ImplicitEdge::interpolate_at(s2, t2, p.y)
    };
}

Cell Cell::subcell(Point const& start, Point const& end) const {
    auto const [start1, start2] = interpolate_at(start);
    auto const [end1, end2] = interpolate_at(end);
    return {start1, start2, end1, end2};
}

/**
 * @brief Computes the integral int_{lo}^{hi} ax + b dx
 * 
 * @param lo 
 * @param hi 
 * @param a 
 * @param b 
 * @return the value of the integral
 */
distance_t integrate_linear_function(distance_t lo, distance_t hi, distance_t a, distance_t b) {
    // std::cout << 0.5*a*pow(hi, 2) + b*hi - (0.5*a*pow(lo, 2) + b*lo) << std::endl;
    return 0.5*a*pow(hi, 2) + b*hi - (0.5*a*pow(lo, 2) + b*lo);
}

distance_t integrate_linear_cost_l1(Cell const& cell,
    Point const& s, Point const& t) {
    auto const [s1, s2] = cell.interpolate_at(s);
    auto const [t1, t2] = cell.interpolate_at(t);

    auto const [dx1, dy1] = s1 - s2;
    auto const [dx2, dy2] = t1 - t2;

    /**
     * 
     * term1 = int_{0}^{1} |s1.x - t*(s1.x - t1.x) - (s2.x - t*(s2.x - t2.x))|dt
     *       = int_{0}^{1} |t*(t1.x-s1.x + s2.x - t2.x) + s1.x - s2.x|dt
     *       
     * 
     * term2 = int_{0}^{1} |s1.y - t*(s1.y - t1.y) - (s2.y - t*(s2.y - t2.y))|dt
     * 
     */

    enum CASE {POS, NEG, POS_NEG, NEG_POS};

    CASE x_term_case;
    CASE y_term_case;

    distance_t x_int;
    distance_t y_int;

    distance_t x_term_at_0 = s1.x - s2.x;
    distance_t x_term_at_1 = t1.x - t2.x;

    distance_t y_term_at_0 = s1.y - s2.y;
    distance_t y_term_at_1 = t1.y - t2.y;

    if ((x_term_at_0 > 0 || approx_zero(x_term_at_0)) && (x_term_at_1 > 0 || approx_zero(x_term_at_1))) {
        x_term_case = POS;
    } else if ((x_term_at_0 < 0 || approx_zero(x_term_at_0)) && (x_term_at_1 < 0 || approx_zero(x_term_at_1))) {
        x_term_case = NEG;
    } else if ((x_term_at_0 > 0 || approx_zero(x_term_at_0)) && (x_term_at_1 < 0 || approx_zero(x_term_at_1))) {
        x_term_case = POS_NEG;
    } else if ((x_term_at_0 < 0 || approx_zero(x_term_at_0)) && (x_term_at_1 > 0 || approx_zero(x_term_at_1))) {
        x_term_case = NEG_POS;
    }

    if ((y_term_at_0 > 0 || approx_zero(y_term_at_0)) && (y_term_at_1 > 0 || approx_zero(y_term_at_1))) {
        y_term_case = POS;
    } else if ((y_term_at_0 < 0 || approx_zero(y_term_at_0)) && (y_term_at_1 < 0 || approx_zero(y_term_at_1))) {
        y_term_case = NEG;
    } else if ((y_term_at_0 > 0 || approx_zero(y_term_at_0)) && (y_term_at_1 < 0 || approx_zero(y_term_at_1))) {
        y_term_case = POS_NEG;
    } else if ((y_term_at_0  < 0|| approx_zero(y_term_at_0)) && (y_term_at_1 > 0 || approx_zero(x_term_at_1))) {
        y_term_case = NEG_POS;
    }

    distance_t xterm_pos_a = t1.x - s1.x + s2.x - t2.x;
    distance_t xterm_pos_b = s1.x - s2.x;
    distance_t yterm_pos_a = t1.y - s1.y + s2.y - t2.y;
    distance_t yterm_pos_b = s1.y - s2.y;

    distance_t xterm_break = !approx_zero(xterm_pos_a) ? -xterm_pos_b / xterm_pos_a : 0;
    distance_t yterm_break =  !approx_zero(yterm_pos_a) ? -yterm_pos_b / yterm_pos_a : 0;

    switch (x_term_case) {
        case POS: {
            x_int = integrate_linear_function(0, 1, xterm_pos_a, xterm_pos_b);
            break;
        }
        case NEG: {
            x_int = integrate_linear_function(0, 1, -xterm_pos_a, -xterm_pos_b);
            break;
        }
        case POS_NEG: {
            distance_t left = integrate_linear_function(0, xterm_break, xterm_pos_a, xterm_pos_b);
            distance_t right = integrate_linear_function(xterm_break, 1, -xterm_pos_a, -xterm_pos_b);
            x_int = left + right;
            break;
        }
        case NEG_POS: {
            distance_t left = integrate_linear_function(0, xterm_break, -xterm_pos_a, -xterm_pos_b);
            distance_t right = integrate_linear_function(xterm_break, 1, xterm_pos_a, xterm_pos_b);
            x_int = left + right;
            break;
        }
    }

    switch (y_term_case) {
        case POS: {
            y_int = integrate_linear_function(0, 1, yterm_pos_a, yterm_pos_b);
            break;
        }
        case NEG: {
            y_int = integrate_linear_function(0, 1, -yterm_pos_a, -yterm_pos_b);
            break;
        }
        case POS_NEG: {
            distance_t left = integrate_linear_function(0, yterm_break, yterm_pos_a, yterm_pos_b);
            distance_t right = integrate_linear_function(yterm_break, 1, -yterm_pos_a, -yterm_pos_b);
            y_int = left + right;
            break;
        }
        case NEG_POS: {
            distance_t left = integrate_linear_function(0, yterm_break, -yterm_pos_a, -yterm_pos_b);
            distance_t right = integrate_linear_function(yterm_break, 1, yterm_pos_a, yterm_pos_b);
            y_int = left + right;
            break;
        }
    }

    return x_int + y_int;
}

// TODO: Move methods out of Cell namespace
distance_t integrate_linear_cost(Cell const& cell,
        Point const& s, Point const& t) {
    auto const [s1, s2] = cell.interpolate_at(s);
    auto const [t1, t2] = cell.interpolate_at(t);

    // Get difference in x- and y-coordinates at the start and end of linear
    // matching (s, t)
    auto const [dx1, dy1] = s1 - s2;
    auto const [dx2, dy2] = t1 - t2;

    distance_t const a = (dx1 - dx2) * (dx1 - dx2) + (dy1 - dy2) * (dy1 - dy2);
    distance_t const b = 2 * (dx1 * dx2 - dx1 * dx1 + dy1 * dy2 - dy1 * dy1);
    distance_t const c = dx1 * dx1 + dy1 * dy1;

    return (a / 3 + b / 2 + c);
}
