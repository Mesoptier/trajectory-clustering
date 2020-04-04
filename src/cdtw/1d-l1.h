#ifndef TRAJECTORY_CLUSTERING_1D_L1_H
#define TRAJECTORY_CLUSTERING_1D_L1_H

#include "cdtw.h"
#include "BivariatePolynomial.h"

//
// 1D + L1 image norm + L1 param norm
//

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
PiecewisePolynomial<2> CDTW<1, Norm::L1, Norm::L1>::base_bottom(const Cell& cell) const {
    // TODO: Clean this up; sx, sy, tx, ty are assumed to be coordinates in a space where (0,0) is the ellipse center
    double sx = cell.s.x - cell.mid.x;
    double sy = cell.s.y - cell.mid.y;
    double tx = cell.t.x - cell.mid.x;

    if (sx >= sy) {
        return PiecewisePolynomial<2>({
            {{sx + cell.mid.x, tx + cell.mid.x}, Polynomial<2>({-(sx * sx) / 2 + sx * sy, -sy, 1. / 2})}
        });
    } else if (sy >= tx) {
        return PiecewisePolynomial<2>({
            {{sx + cell.mid.x, tx + cell.mid.x}, Polynomial<2>({(sx * sx) / 2 - sx * sy, sy, -1. / 2})}
        });
    } else { // sx < sy < tx
        return PiecewisePolynomial<2>({
            {{sx + cell.mid.x, sy + cell.mid.x}, Polynomial<2>({(sx * sx) / 2 - sx * sy, sy, -1. / 2})},
            {{sy + cell.mid.x, tx + cell.mid.x}, Polynomial<2>({(sx * sx) / 2 - sx * sy + sy * sy, -sy, 1. / 2})}
        });
    }
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
CDTW<1, Norm::L1, Norm::L1>::bottom_to_right_costs(const Cell& cell) const {
    // Coordinates of the cell in a system where cell.mid lies on (0, 0).
    // This makes it much easier to formulate the bivariate polynomials and constraints, but does require us to
    // translate the cell back such that (s.x, s.y) = (0, 0) using .translate_xy(-sx, -sy).
    double sx = cell.s.x - cell.mid.x;
    double sy = cell.s.y - cell.mid.y;
    double tx = cell.t.x - cell.mid.x;
    double ty = cell.t.y - cell.mid.y;

    std::vector<ConstrainedBivariatePolynomial<2>> costs;
    costs.reserve(5);

    // 1.
    // Not via valley
    // A below valley
    // B below valley
    costs.push_back(ConstrainedBivariatePolynomial<2>{
        BivariatePolynomial<2>({{
            {{(sy * sy) / 2 + (tx * tx) / 2, -tx, -1./2}},
            {{-sy, 2, 0}},
            {{-1./2, 0, 0}},
        }}),
        // B in cell + B below valley
        Interval{sy, std::clamp(tx, sy, ty)},
        {{
            Polynomial<1>({sx, 0}), // A in cell
            Polynomial<1>({sy, 0}), // A below valley
            Polynomial<1>({0, 1}), // Not via valley
        }},
        {{
             Polynomial<1>({tx, 0}), // A in cell
        }}
    }.translate_xy(-sx, -sy));

    // 2.
    // Via valley
    // A below valley
    // B below valley
    costs.push_back(ConstrainedBivariatePolynomial<2>{
        BivariatePolynomial<2>({{
            {{(sy * sy) / 2 + (tx * tx) / 2, -tx, 1./2}},
            {{-sy, 0, 0}},
            {{1./2, 0, 0}},
        }}),
        // B in cell + B below valley
        Interval{sy, std::clamp(tx, sy, ty)},
        {{
            Polynomial<1>({sx, 0}), // A in cell
            Polynomial<1>({sy, 0}), // A below valley
        }},
        {{
             Polynomial<1>({tx, 0}), // A in cell
             Polynomial<1>({0, 1}), // Via valley
        }}
    }.translate_xy(-sx, -sy));

    // 3.
    // Via valley
    // A above valley
    // B below valley
    costs.push_back(ConstrainedBivariatePolynomial<2>{
        BivariatePolynomial<2>({{
            {{(sy * sy) / 2 + (tx * tx) / 2, -tx, 1./2}},
            {{-sy, 0, 0}},
            {{1./2, 0, 0}},
        }}),
        // B in cell + B below valley
        Interval{sy, std::clamp(tx, sy, ty)},
        {{
            Polynomial<1>({sx, 0}), // A in cell
        }},
        {{
             Polynomial<1>({tx, 0}), // A in cell
             Polynomial<1>({sy, 0}), // A above valley
        }}
    }.translate_xy(-sx, -sy));

    // 4.
    // Via valley
    // A below valley
    // B above valley
    costs.push_back(ConstrainedBivariatePolynomial<2>{
        BivariatePolynomial<2>({{
            {{(sy * sy) / 2 + (tx * tx) / 2, -tx, 1./2}},
            {{-sy, 0, 0}},
            {{1./2, 0, 0}},
        }}),
        // B in cell + B above valley
        Interval{std::clamp(tx, sy, ty), ty},
        {{
            Polynomial<1>({sx, 0}), // A in cell
            Polynomial<1>({sy, 0}), // A below valley
        }},
        {{
             Polynomial<1>({tx, 0}), // A in cell
        }}
    }.translate_xy(-sx, -sy));

    // 5. / 6.
    // A above valley
    // B above valley
    if (tx > sy) { // Via valley
        costs.push_back(ConstrainedBivariatePolynomial<2>{
            BivariatePolynomial<2>({{
                {{(sy * sy) / 2 + (tx * tx) / 2, -tx, 1./2}},
                {{-sy, 0, 0}},
                {{1./2, 0, 0}},
            }}),
            // B in cell + B above valley
            Interval{std::clamp(tx, sy, ty), ty},
            {{
                Polynomial<1>({sx, 0}), // A in cell
            }},
            {{
                 Polynomial<1>({tx, 0}), // A in cell
                 Polynomial<1>({sy, 0}), // A above valley
            }}
        }.translate_xy(-sx, -sy));
    } else { // Not via valley
        costs.push_back(ConstrainedBivariatePolynomial<2>{
            BivariatePolynomial<2>({{
                {{2 * sy * tx - (sy * sy) / 2 - (tx * tx) / 2, -tx, 1./2}},
                {{-sy, 0, 0}},
                {{1./2, 0, 0}},
            }}),
            // B in cell + B above valley
            Interval{std::clamp(tx, sy, ty), ty},
            {{
                Polynomial<1>({sx, 0}), // A in cell
            }},
            {{
                 Polynomial<1>({tx, 0}), // A in cell
                 Polynomial<1>({sy, 0}), // A above valley
            }}
        }.translate_xy(-sx, -sy));
    }

    return costs;
}

PiecewisePolynomial<2> bottom_to_top_6(double sx, double sy, double tx, double ty) {
    // Not via valley
    // A below valley
    // B below valley

    BivariatePolynomial<2> h({{
        {{(sy * sy) / 2 - (ty * ty) / 2, -ty, 1./2}},
        {{-sy + 2 * ty, 0, 0}},
        {{-1./2, 0, 0}},
    }});

    // B in cell + B below valley
    Interval b_interval{std::clamp(ty, sx, tx), tx};

    std::vector<Polynomial<1>> left_constraints;
    std::vector<Polynomial<1>> right_constraints;

    // A in cell
    left_constraints.push_back(Polynomial<1>({sx, 0}));
    right_constraints.push_back(Polynomial<1>({tx, 0}));

    // B to the right of A
    right_constraints.push_back(Polynomial<1>({0, 1}));

    // A below valley
    left_constraints.push_back(Polynomial<1>({sy, 0}));

    // Not via valley
    left_constraints.push_back(Polynomial<1>({ty, 0}));

    return find_minimum(h, b_interval, left_constraints, right_constraints);
}

PiecewisePolynomial<2> bottom_to_top_5(double sx, double sy, double tx, double ty) {
    // Via valley
    // A below valley
    // B below valley

    BivariatePolynomial<2> h({{
        {{(sy * sy) / 2 + (ty * ty) / 2, -ty, 1./2}},
        {{-sy, 0, 0}},
        {{1./2, 0, 0}},
    }});

    // B in cell + B below valley
    Interval b_interval{std::clamp(ty, sx, tx), tx};

    std::vector<Polynomial<1>> left_constraints;
    std::vector<Polynomial<1>> right_constraints;

    // A in cell
    left_constraints.push_back(Polynomial<1>({sx, 0}));
    right_constraints.push_back(Polynomial<1>({tx, 0}));

    // B to the right of A
    right_constraints.push_back(Polynomial<1>({0, 1}));

    // A below valley
    left_constraints.push_back(Polynomial<1>({sy, 0}));

    // Via valley
    left_constraints.push_back(Polynomial<1>({ty, 0}));

    return find_minimum(h, b_interval, left_constraints, right_constraints);
}

PiecewisePolynomial<2> bottom_to_top_4(double sx, double sy, double tx, double ty) {
    // Via valley
    // A below valley
    // B above valley

    BivariatePolynomial<2> h({{
        {{(sy * sy) / 2 + (ty * ty) / 2, -ty, 1./2}},
        {{-sy, 0, 0}},
        {{1./2, 0, 0}},
    }});

    // B in cell + B above valley
    Interval b_interval{sx, std::clamp(ty, sx, tx)};

    std::vector<Polynomial<1>> left_constraints;
    std::vector<Polynomial<1>> right_constraints;

    // A in cell
    left_constraints.push_back(Polynomial<1>({sx, 0}));
    right_constraints.push_back(Polynomial<1>({tx, 0}));

    // B to the right of A
    right_constraints.push_back(Polynomial<1>({0, 1}));

    // A below valley
    left_constraints.push_back(Polynomial<1>({sy, 0}));

    // Via valley
    right_constraints.push_back(Polynomial<1>({ty, 0}));

    return find_minimum(h, b_interval, left_constraints, right_constraints);
}

PiecewisePolynomial<2> bottom_to_top_3(double sx, double sy, double tx, double ty) {
    // Via valley
    // A above valley
    // B below valley

    BivariatePolynomial<2> h({{
        {{(sy * sy) / 2 + (ty * ty) / 2, -ty, 1./2}},
        {{-sy, 0, 0}},
        {{1./2, 0, 0}},
    }});

    // B in cell + B below valley
    Interval b_interval{std::clamp(ty, sx, tx), tx};

    std::vector<Polynomial<1>> left_constraints;
    std::vector<Polynomial<1>> right_constraints;

    // A in cell
    left_constraints.push_back(Polynomial<1>({sx, 0}));
    right_constraints.push_back(Polynomial<1>({tx, 0}));

    // B to the right of A
    right_constraints.push_back(Polynomial<1>({0, 1}));

    // A above valley
    right_constraints.push_back(Polynomial<1>({sy, 0}));

    return find_minimum(h, b_interval, left_constraints, right_constraints);
}

PiecewisePolynomial<2> bottom_to_top_2(double sx, double sy, double tx, double ty) {
    // Via valley
    // A above valley
    // B above valley

    BivariatePolynomial<2> h({{
        {{(sy * sy) / 2 + (ty * ty) / 2, -ty, 1./2}},
        {{-sy, 0, 0}},
        {{1./2, 0, 0}},
    }});

    // B in cell + B above valley + Via valley
    Interval b_interval{std::clamp(sy, sx, tx), std::clamp(ty, sx, tx)};

    std::vector<Polynomial<1>> left_constraints;
    std::vector<Polynomial<1>> right_constraints;

    // A in cell
    left_constraints.push_back(Polynomial<1>({sx, 0}));
    right_constraints.push_back(Polynomial<1>({tx, 0}));

    // B to the right of A
    right_constraints.push_back(Polynomial<1>({0, 1}));

    // A above valley
    right_constraints.push_back(Polynomial<1>({sy, 0}));

    return find_minimum(h, b_interval, left_constraints, right_constraints);
}

PiecewisePolynomial<2> bottom_to_top_1(double sx, double sy, double tx, double ty) {
    // Not via valley
    // A above valley
    // B above valley

    BivariatePolynomial<2> h({{
        {{(ty * ty) / 2 - (sy * sy) / 2, 2 * sy - ty, -1./2}},
        {{-sy, 0, 0}},
        {{1./2, 0, 0}},
    }});

    // B in cell + B above valley + Via valley
    Interval b_interval{sx, std::clamp(sy, sx, tx)};

    std::vector<Polynomial<1>> left_constraints;
    std::vector<Polynomial<1>> right_constraints;

    // A in cell
    left_constraints.push_back(Polynomial<1>({sx, 0}));
    right_constraints.push_back(Polynomial<1>({tx, 0}));

    // B to the right of A
    right_constraints.push_back(Polynomial<1>({0, 1}));

    // A above valley
    right_constraints.push_back(Polynomial<1>({sy, 0}));

    return find_minimum(h, b_interval, left_constraints, right_constraints);
}

PiecewisePolynomial<2> bottom_to_top(double sx, double sy, double tx, double ty) {
    auto result = bottom_to_top_1(sx, sy, tx, ty);
    fast_lower_envelope(result, bottom_to_top_2(sx, sy, tx, ty));
    fast_lower_envelope(result, bottom_to_top_3(sx, sy, tx, ty));
    fast_lower_envelope(result, bottom_to_top_4(sx, sy, tx, ty));
    fast_lower_envelope(result, bottom_to_top_5(sx, sy, tx, ty));
    fast_lower_envelope(result, bottom_to_top_6(sx, sy, tx, ty));
    return result;
}

#endif //TRAJECTORY_CLUSTERING_1D_L1_H
