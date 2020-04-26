#ifndef TRAJECTORY_CLUSTERING_1D_L2SQUARED_L1_H
#define TRAJECTORY_CLUSTERING_1D_L2SQUARED_L1_H

#include "cdtw.h"
#include "BivariatePolynomial.h"

//
// 1D + L2Squared image norm + L1 param norm
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
PiecewisePolynomial<3> CDTW<1, Norm::L2Squared, Norm::L1>::base_bottom(const Cell& cell) const {
    // TODO: Clean this up; sx, sy, tx, ty are assumed to be coordinates in a space where (0,0) is the ellipse center
    const double sx = cell.s.x - cell.mid.x;
    const double sy = cell.s.y - cell.mid.y;
    const double tx = cell.t.x - cell.mid.x;

    bool same_direction = (cell.s1.x <= cell.t1.x) == (cell.s2.x <= cell.t2.x);

    if (!same_direction) { // Downwards valley
        return PiecewisePolynomial<3>({
            {{0, tx-sx}, Polynomial<3>({0, (sx+sy)*(sx+sy), sx+sy, 1./3})},
        });
    }

    // Upwards valley
    return PiecewisePolynomial<3>({
        {{0, tx-sx}, Polynomial<3>({0, (sx-sy)*(sx-sy), sx-sy, 1./3})},
    });
}

/**
 * Compute the piecewise bivariate polynomial representing the cost of the optimal path between a point A on the bottom
 * boundary and a point B on the right boundary.
 *
 * @param cell
 * @return A piecewise bivariate polynomial over the domain 0 <= A <= cell.width and 0 <=B <= cell.height.
 */
template<>
std::vector<ConstrainedBivariatePolynomial<3>>
CDTW<1, Norm::L2Squared, Norm::L1>::bottom_to_right_costs(const Cell& cell) const {
    // Coordinates of the cell in a system where cell.mid lies on (0, 0).
    // This makes it much easier to formulate the bivariate polynomials and constraints, but does require us to
    // translate the cell back such that (s.x, s.y) = (0, 0) using .translate_xy(-sx, -sy).
    const double sx = cell.s.x - cell.mid.x;
    const double sy = cell.s.y - cell.mid.y;
    const double tx = cell.t.x - cell.mid.x;
    const double ty = cell.t.y - cell.mid.y;

    std::vector<ConstrainedBivariatePolynomial<3>> costs;

    bool same_direction = (cell.s1.x <= cell.t1.x) == (cell.s2.x <= cell.t2.x);

    if (!same_direction) { // Downwards valley
        costs.reserve(1);

        // 1.
        costs.push_back(ConstrainedBivariatePolynomial<3>{
            BivariatePolynomial<3>({{
                {{(sx - tx) * (sx*sx + 3*sy*sy + 3*sy*tx + tx*tx + 3*sy*sx + tx*sx) / -3, (sy+tx)*(sy+tx), sy + tx, 1./3}},
                {{-(sx+sy) * (sx+sy), 0, 0, 0}},
                {{-sx - sy, 0, 0, 0}},
                {{-1./3, 0, 0, 0}},
            }}),
            // B in cell
            Interval{0, ty-sy},
            {{
                Polynomial<1>({0, 0}), // A in cell
            }},
            {{
                Polynomial<1>({tx-sx, 0}), // A in cell
            }},
        });

        return costs;
    }

    // Upwards valley

    costs.reserve(5);

    // 1.
    // Not via valley
    // A below valley
    // B below valley
    costs.push_back(ConstrainedBivariatePolynomial<3>{
        BivariatePolynomial<3>({{
            {{
                (sx - tx) * (sx*sx + 3*sy*sy - 3*sy*tx + tx*tx - 3*sy*sx + tx*sx) / -3,
                2*sx*sx - 4*sx*sy + sy*sy + 2*sy*tx - tx*tx,
                -2*sx + sy + tx,
                1./3,
            }},
            {{-(sx-sy) * (sx-sy), 4 * (sx-sy), -2, 0}},
            {{-sx + sy, 2, 0, 0}},
            {{-1./3, 0, 0, 0}},
        }}),
        // B in cell + B below valley
        Interval{0, std::clamp(sx-sy, 0., ty-sy)},
        {{
            Polynomial<1>({0, 0}), // A in cell
            Polynomial<1>({sy-sx, 0}), // A below valley
            Polynomial<1>({sy-sx, 1}), // Not via valley
        }},
        {{
            Polynomial<1>({tx-sx, 0}), // A in cell
        }}
    });

    // 2.
    // Via valley
    // A below valley
    // B below valley
    costs.push_back(ConstrainedBivariatePolynomial<3>{
        BivariatePolynomial<3>({{
            {{
                (sx - 2*sy + tx) * (sx*sx + sy*sy - sy*tx + tx*tx - sx*sy - sx*tx) / 3,
                -(sy-tx) * (sy-tx),
                -sy + tx,
                -1./3,
            }},
            {{(sx-sy) * (sx-sy), 0, 0, 0}},
            {{sx - sy, 0, 0, 0}},
            {{1./3, 0, 0, 0}},
        }}),
        // B in cell + B below valley
        Interval{0, std::clamp(sx-sy, 0., ty-sy)},
        {{
            Polynomial<1>({0, 0}), // A in cell
            Polynomial<1>({sy-sx, 0}), // A below valley
        }},
        {{
            Polynomial<1>({tx-sx, 0}), // A in cell
            Polynomial<1>({sy-sx, 1}), // Via valley
        }}
    });

    // 3.
    // Via valley
    // A above valley
    // B below valley
    costs.push_back(ConstrainedBivariatePolynomial<3>{
        BivariatePolynomial<3>({{
            {{
                (sx - tx) * (sx*sx + 3*sy*sy - 3*sy*tx + tx*tx - 3*sy*sx + tx*sx) / -3,
                -(sy-tx) * (sy-tx),
                -sy + tx,
                -1./3,
            }},
            {{-(sx-sy) * (sx-sy), 0, 0, 0}},
            {{-sx + sy, 0, 0, 0}},
            {{-1./3, 0, 0, 0}},
        }}),
        // B in cell + B below valley
        Interval{0, std::clamp(sx-sy, 0., ty-sy)},
        {{
            Polynomial<1>({0, 0}), // A in cell
        }},
        {{
            Polynomial<1>({tx-sx, 0}), // A in cell
            Polynomial<1>({sy-sx, 0}), // A above valley
        }}
    });

    // 4.
    // Via valley
    // A below valley
    // B above valley
    costs.push_back(ConstrainedBivariatePolynomial<3>{
        BivariatePolynomial<3>({{
            {{
                (sx - tx) * (sx*sx + 3*sy*sy - 3*sy*tx + tx*tx - 3*sy*sx + tx*sx) / 3,
                (sy-tx) * (sy-tx),
                sy - tx,
                1./3
            }},
            {{(sx-sy) * (sx-sy), 0, 0, 0}},
            {{sx - sy, 0, 0, 0}},
            {{1./3, 0, 0, 0}},
        }}),
        // B in cell + B above valley
        Interval{std::clamp(sx-sy, 0., ty-sy), ty-sy},
        {{
            Polynomial<1>({0, 0}), // A in cell
            Polynomial<1>({sy-sx, 0}), // A below valley
        }},
        {{
            Polynomial<1>({tx-sx, 0}), // A in cell
        }}
    });

    // 5. / 6.
    // A above valley
    // B above valley
    if (tx > sy) { // Via valley
        costs.push_back(ConstrainedBivariatePolynomial<3>{
            BivariatePolynomial<3>({{
                {{
                    (sx - 2*sy + tx) * (sx*sx + sy*sy - sy*tx + tx*tx - sy*sx + tx*sx) / -3,
                    (sy-tx) * (sy-tx),
                    sy - tx,
                    1./3,
                }},
                {{-(sx-sy) * (sx-sy), 0, 0, 0}},
                {{-sx + sy, 0, 0, 0}},
                {{-1./3, 0, 0, 0}},
            }}),
            // B in cell + B above valley
            Interval{std::clamp(sx-sy, 0., ty-sy), ty-sy},
            {{
                Polynomial<1>({0, 0}), // A in cell
            }},
            {{
                Polynomial<1>({tx-sx, 0}), // A in cell
                Polynomial<1>({sy-sx, 0}), // A below valley
            }}
        });
    } else { // Not via valley
        costs.push_back(ConstrainedBivariatePolynomial<3>{
            BivariatePolynomial<3>({{
                {{
                    ((tx-ty)*(tx-ty)*(tx-ty) - (sx-sy)*(sx-sy)*(sx-sy)) / 3,
                    (sy-tx)*(sy-tx),
                    sy - tx,
                    1./3,
                }},
                {{-(sx-sy)*(sx-sy), 0, 0, 0}},
                {{-sx+sy, 0, 0, 0}},
                {{-1./3, 0, 0, 0}},
            }}),
            // B in cell + B above valley
            Interval{std::clamp(sx-sy, 0., ty-sy), ty-sy},
            {{
                Polynomial<1>({0, 0}), // A in cell
            }},
            {{
                Polynomial<1>({tx-sx, 0}), // A in cell
                Polynomial<1>({sy-sx, 0}), // A below valley
            }}
        });
    }

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
std::vector<ConstrainedBivariatePolynomial<3>>
CDTW<1, Norm::L2Squared, Norm::L1>::bottom_to_top_costs(const Cell& cell) const {
    // Coordinates of the cell in a system where cell.mid lies on (0, 0).
    // This makes it much easier to formulate the bivariate polynomials and constraints, but does require us to
    // translate the cell back such that (s.x, s.y) = (0, 0) using .translate_xy(-sx, -sy).
    const double sx = cell.s.x - cell.mid.x;
    const double sy = cell.s.y - cell.mid.y;
    const double tx = cell.t.x - cell.mid.x;
    const double ty = cell.t.y - cell.mid.y;

    std::vector<ConstrainedBivariatePolynomial<3>> costs;

    bool same_direction = (cell.s1.x <= cell.t1.x) == (cell.s2.x <= cell.t2.x);

    if (!same_direction) { // Downwards valley
        costs.reserve(1);

        // 1.
        costs.push_back(ConstrainedBivariatePolynomial<3>{
            BivariatePolynomial<3>({{
                {{(sy - ty) * (3 * sx*sx + sy*sy + sy*ty + ty*ty - 3*sx*(sy+ty)) / -3, sx*sx + 2*sy*sy - ty*ty + 2*sx*(-2*sy+ty), sx -2*sy + ty, 1./3}},
                {{-(sx-sy) * (sx-sy), 0, 0, 0}},
                {{-sx + sy, 0, 0, 0}},
                {{-1./3, 0, 0, 0}},
            }}),
            // B in cell
            Interval{0, tx-sx},
            {{
                Polynomial<1>({0, 0}), // A in cell
            }},
            {{
                Polynomial<1>({tx-sx, 0}), // A in cell
                Polynomial<1>({0, 1}), // A <= B
            }},
        });

        return costs;
    }

    // Upwards valley

    costs.reserve(6);

    // 1.
    // Not via valley
    // A above valley
    // B above valley
    costs.push_back(ConstrainedBivariatePolynomial<3>{
        BivariatePolynomial<3>({{
            {{
                (sy - ty) * (3*sx*sx + sy*sy + sy*ty + ty*ty - 3*sx*sy - 3*sx*ty) / -3,
                sx*sx + 2*sy*sy - ty*ty + 2*sx*(ty - 2*sy),
                sx - 2*sy + ty,
                1./3,
            }},
            {{-(sx-sy) * (sx-sy), 0, 0, 0}},
            {{-sx + sy, 0, 0, 0}},
            {{-1./3, 0, 0, 0}},
        }}),
        // B in cell + B above valley + Not via valley
        Interval{0, std::clamp(sy-sx, 0., tx-sx)},
        {{
            Polynomial<1>({0, 0}), // A in cell
        }},
        {{
            Polynomial<1>({tx-sx, 0}), // A in cell
            Polynomial<1>({0, 1}), // A <= B
            Polynomial<1>({sy-sx, 0}), // A above valley
        }}
    });

    // 2.
    // Via valley
    // A above valley
    // B above valley
    costs.push_back(ConstrainedBivariatePolynomial<3>{
        BivariatePolynomial<3>({{
            {{
                (-2*sx*sx*sx + sy*sy*sy + ty*ty*ty + 3*sx*sx*(sy+ty) - 3*sx*(sy*sy + ty*ty)) / 3,
                -(sx-ty) * (sx-ty),
                -sx + ty,
                -1./3,
            }},
            {{-(sx-sy) * (sx-sy), 0, 0, 0}},
            {{-sx + sy, 0, 0, 0}},
            {{-1./3, 0, 0, 0}},
        }}),
        // B in cell + B above valley + Via valley
        Interval{std::clamp(sy-sx, 0., tx-sx), std::clamp(ty-sx, 0., tx-sx)},
        {{
            Polynomial<1>({0, 0}), // A in cell
        }},
        {{
            Polynomial<1>({tx-sx, 0}), // A in cell
            Polynomial<1>({0, 1}), // A <= B
            Polynomial<1>({sy-sx, 0}), // A above valley
        }}
    });

    // 3.
    // Via valley
    // A below valley
    // B above valley
    costs.push_back(ConstrainedBivariatePolynomial<3>{
        BivariatePolynomial<3>({{
            {{
                (sy - ty) * (3*sx*sx + sy*sy + sy*ty + ty*ty - 3*sx*sy - 3*sx*ty) / -3,
                -(sx-ty) * (sx-ty),
                -sx + ty,
                -1./3,
            }},
            {{(sx-sy)*(sx-sy), 0, 0, 0}},
            {{sx - sy, 0, 0, 0}},
            {{1./3, 0, 0, 0}},
        }}),
        // B in cell + B above valley
        Interval{0, std::clamp(ty-sx, 0., tx-sx)},
        {{
            Polynomial<1>({0, 0}), // A in cell
            Polynomial<1>({sy-sx, 0}), // A below valley
        }},
        {{
            Polynomial<1>({tx-sx, 0}), // A in cell
            Polynomial<1>({0, 1}), // A <= B
        }}
    });

    // 4.
    // Via valley
    // A above valley
    // B below valley
    costs.push_back(ConstrainedBivariatePolynomial<3>{
        BivariatePolynomial<3>({{
            {{
                (sy - ty) * (3*sx*sx + sy*sy + sy*ty + ty*ty - 3*sx*sy - 3*sx*ty) / 3,
                (sx-ty) * (sx-ty),
                sx - ty,
                1./3,
            }},
            {{-(sx-sy) * (sx-sy), 0, 0, 0}},
            {{-sx + sy, 0, 0, 0}},
            {{-1./3, 0, 0, 0}},
        }}),
        // B in cell + B below valley
        Interval{std::clamp(ty-sx, 0., tx-sx), tx-sx},
        {{
            Polynomial<1>({0, 0}), // A in cell
        }},
        {{
            Polynomial<1>({tx-sx, 0}), // A in cell
            Polynomial<1>({0, 1}), // A <= B
            Polynomial<1>({sy-sx, 0}), // A above valley
        }}
    });

    // 5.
    // Via valley
    // A below valley
    // B below valley
    costs.push_back(ConstrainedBivariatePolynomial<3>{
        BivariatePolynomial<3>({{
            {{
                (2*sx - sy - ty) * (sx*sx + sy*sy - sy*ty + ty*ty - sx*sy - sx*ty) / 3,
                (sx-ty) * (sx-ty),
                sx - ty,
                1./3
            }},
            {{(sx-sy) * (sx-sy), 0, 0, 0}},
            {{sx - sy, 0, 0, 0}},
            {{1./3, 0, 0, 0}},
        }}),
        // B in cell + B below valley
        Interval{std::clamp(ty-sx, 0., tx-sx), tx-sx},
        {{
            Polynomial<1>({0, 0}), // A in cell
            Polynomial<1>({sy-sx, 0}), // A below valley
        }},
        {{
            Polynomial<1>({tx-sx, 0}), // A in cell
            Polynomial<1>({0, 1}), // A <= B
            Polynomial<1>({ty-sx, 0}), // Via valley
        }}
    });

    // 6.
    // Not via valley
    // A below valley
    // B below valley
    costs.push_back(ConstrainedBivariatePolynomial<3>{
        BivariatePolynomial<3>({{
            {{
                (sy - ty) * (3*sx*sx + sy*sy + sy*ty + ty*ty - 3*sx*sy - 3*sx*ty) / 3,
                (sx-ty) * (sx-ty),
                sx - ty,
                1./3,
            }},
            {{-sx*sx - 2*sx*sy + sy*sy + 4*sx*ty - 2*ty*ty, 0, 0, 0}},
            {{-sx - sy + 2*ty, 0, 0, 0}},
            {{-1./3, 0, 0, 0}},
        }}),
        // B in cell + B below valley
        Interval{std::clamp(ty-sx, 0., tx-sx), tx-sx},
        {{
            Polynomial<1>({0, 0}), // A in cell
            Polynomial<1>({ty-sx, 0}), // A below valley + Not via valley
        }},
        {{
            Polynomial<1>({tx-sx, 0}), // A in cell
            Polynomial<1>({0, 1}), // A <= B
        }}
    });

    return costs;
}

#endif //TRAJECTORY_CLUSTERING_1D_L2SQUARED_L1_H
