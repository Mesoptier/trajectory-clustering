#ifndef TRAJECTORY_CLUSTERING_1D_L1_L1_H
#define TRAJECTORY_CLUSTERING_1D_L1_L1_H

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

    bool same_direction = (cell.s1.x <= cell.t1.x) == (cell.s2.x <= cell.t2.x);

    if (!same_direction) { // Downwards valley
        if (sx >= -sy) {
            return PiecewisePolynomial<2>({
                {{0, tx-sx}, Polynomial<2>({0, sx + sy, 1./2})},
            });
        } else if (-sy >= tx) {
            return PiecewisePolynomial<2>({
                {{0, tx-sx}, Polynomial<2>({0, -sx-sy, -1./2})},
            });
        } else { // sx < -sy < tx
            return PiecewisePolynomial<2>({
                {{0, -sy-sx}, Polynomial<2>({0, -sx-sy, -1./2})},
                {{-sy-sx, tx-sx}, Polynomial<2>({(sx+sy) * (sx+sy), sx+sy, 1./2})},
            });
        }
    }

    // Upwards valley

    if (sy <= sx) {
        return PiecewisePolynomial<2>({
            {{sx, tx}, Polynomial<2>({-(sx * sx) / 2 + sx * sy, -sy, 1. / 2})}
        }).translate(-sx);
    } else if (sy >= tx) {
        return PiecewisePolynomial<2>({
            {{sx, tx}, Polynomial<2>({(sx * sx) / 2 - sx * sy, sy, -1. / 2})}
        }).translate(-sx);
    } else { // sx < sy < tx
        return PiecewisePolynomial<2>({
            {{sx, sy}, Polynomial<2>({(sx * sx) / 2 - sx * sy, sy, -1. / 2})},
            {{sy, tx}, Polynomial<2>({(sx * sx) / 2 - sx * sy + sy * sy, -sy, 1. / 2})}
        }).translate(-sx);
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

    bool same_direction = (cell.s1.x <= cell.t1.x) == (cell.s2.x <= cell.t2.x);

    if (!same_direction) { // Downwards valley
        costs.reserve(3);

        // 1.
        // A below valley
        // Corner below valley
        // B below valley
        costs.push_back(ConstrainedBivariatePolynomial<2>{
            BivariatePolynomial<2>({{
                {{(sx * sx - tx * tx) / 2 + (sx - tx) * sy, -sy - tx, -1./2}},
                {{sx + sy, 0, 0}},
                {{1./2, 0, 0}},
            }}),
            // B in cell + B below valley
            Interval{0, std::clamp(-tx-sy, 0., ty-sy)},
            {{
                Polynomial<1>({0, 0}), // A in cell
            }},
            {{
                 Polynomial<1>({tx-sx, 0}), // A in cell
                 Polynomial<1>({-sy-sx, 0}), // A below valley
            }},
        });

        // 2.
        // A below valley
        // Corner below/above valley
        // B above valley
        costs.push_back(ConstrainedBivariatePolynomial<2>{
            BivariatePolynomial<2>({{
                {{(sx * sx + tx * tx) / 2 + (sx + tx + sy) * sy, sy + tx, 1./2}},
                {{sx + sy, 0, 0}},
                {{1./2, 0, 0}},
            }}),
            // B in cell + B above valley
            Interval{std::clamp(-tx-sy, 0., ty-sy), ty-sy},
            {{
                 Polynomial<1>({0, 0}), // A in cell
             }},
            {{
                 Polynomial<1>({tx-sx, 0}), // A in cell
                 Polynomial<1>({-sy-sx, 0}), // A below valley
             }},
        });

        // 3.
        // A above valley
        // Corner above valley
        // B above valley
        costs.push_back(ConstrainedBivariatePolynomial<2>{
            BivariatePolynomial<2>({{
                {{(tx * tx + sx * sx) / 2 + (tx - sx) * sy, sy + tx, 1./2}},
                {{-sx - sy, 0, 0}},
                {{-1./2, 0, 0}},
            }}),
            // B in cell + B above valley
            Interval{std::clamp(-tx-sy, 0., ty-sy), ty-sy},
            {{
                Polynomial<1>({0, 0}), // A in cell
                Polynomial<1>({-sy-sx, 0}), // A above valley
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

/**
 * Compute the piecewise bivariate polynomial representing the cost of the optimal path between a point A on the bottom
 * boundary and a point B on the top boundary.
 *
 * @param cell
 * @return A piecewise bivariate polynomial over the domain 0 <= a <= b <= cell.width.
 */
template<>
std::vector<ConstrainedBivariatePolynomial<2>>
CDTW<1, Norm::L1, Norm::L1>::bottom_to_top_costs(const Cell& cell) const {
    // Coordinates of the cell in a system where cell.mid lies on (0, 0).
    // This makes it much easier to formulate the bivariate polynomials and constraints, but does require us to
    // translate the cell back such that (s.x, s.y) = (0, 0) using .translate_xy(-sx, -sy).
    double sx = cell.s.x - cell.mid.x;
    double sy = cell.s.y - cell.mid.y;
    double tx = cell.t.x - cell.mid.x;
    double ty = cell.t.y - cell.mid.y;

    std::vector<ConstrainedBivariatePolynomial<2>> costs;

    bool same_direction = (cell.s1.x <= cell.t1.x) == (cell.s2.x <= cell.t2.x);

    if (!same_direction) { // Downwards valley
        costs.reserve(3);

        // 1.
        // A below valley
        // Corner below valley
        // B below valley
        costs.push_back(ConstrainedBivariatePolynomial<2>{
            BivariatePolynomial<2>({{
                {{(sy*sy - ty*ty) / 2 + (sy-ty) * sx, -sx-ty, -1./2}},
                {{sx + sy, 0, 0}},
                {{1./2, 0, 0}},
            }}),
            // B in cell + B below valley
            Interval{0, std::clamp(-ty-sx, 0., tx-sx)},
            {{
                Polynomial<1>({0, 0}), // A in cell
            }},
            {{
                 Polynomial<1>({tx-sx, 0}), // A in cell
                 Polynomial<1>({-sy-sx, 0}), // A below valley
                 Polynomial<1>({0, 1}), // A <= B
            }},
        });

        // 2.
        // A below valley
        // Corner below/above valley
        // B above valley
        costs.push_back(ConstrainedBivariatePolynomial<2>{
            BivariatePolynomial<2>({{
                {{(sy*sy + ty*ty) / 2 + (sx + sy + ty) * sx, sx + ty, 1./2}},
                {{sx + sy, 0, 0}},
                {{1./2, 0, 0}},
            }}),
            // B in cell + B above valley
            Interval{std::clamp(-ty-sx, 0., tx-sx), tx-sx},
            {{
                Polynomial<1>({0, 0}), // A in cell
            }},
            {{
                 Polynomial<1>({tx-sx, 0}), // A in cell
                 Polynomial<1>({-sy-sx, 0}), // A below valley
                 Polynomial<1>({0, 1}), // A <= B
            }},
        });

        // 3.
        // A above valley
        // Corner above valley
        // B above valley
        costs.push_back(ConstrainedBivariatePolynomial<2>{
            BivariatePolynomial<2>({{
                {{(ty*ty - sy*sy) / 2 + (ty - sy) * sx, sx+ty, 1./2}},
                {{-sx-sy, 0, 0}},
                {{-1./2, 0, 0}},
            }}),
            // B in cell + B above valley
            Interval{std::clamp(-ty-sx, 0., tx-sx), tx-sx},
            {{
                Polynomial<1>({0, 0}), // A in cell
                Polynomial<1>({-sy-sx, 0}), // A above valley
            }},
            {{
                 Polynomial<1>({tx-sx, 0}), // A in cell
                 Polynomial<1>({0, 1}), // A <= B
            }},
        });

        return costs;
    }

    // Upwards valley

    costs.reserve(5);

    // 1.
    // Not via valley
    // A above valley
    // B above valley
    costs.push_back(ConstrainedBivariatePolynomial<2>{
        BivariatePolynomial<2>({{
            {{(ty * ty) / 2 - (sy * sy) / 2, 2 * sy - ty, -1./2}},
            {{-sy, 0, 0}},
            {{1./2, 0, 0}},
        }}),
        // B in cell + B above valley + Not via valley
        {sx, std::clamp(sy, sx, tx)},
        {{
             Polynomial<1>({sx, 0}), // A in cell
        }},
        {{
             Polynomial<1>({tx, 0}), // A in cell
             Polynomial<1>({0, 1}), // A <= B
             Polynomial<1>({sy, 0}), // A above valley
        }}
    }.translate_xy(-sx, -sx));

    // 2.
    // Via valley
    // A above valley
    // B above valley
    costs.push_back(ConstrainedBivariatePolynomial<2>{
        BivariatePolynomial<2>({{
            {{(sy * sy) / 2 + (ty * ty) / 2, -ty, 1./2}},
            {{-sy, 0, 0}},
            {{1./2, 0, 0}},
        }}),
        // B in cell + B above valley + Via valley
        {std::clamp(sy, sx, tx), std::clamp(ty, sx, tx)},
        {{
             Polynomial<1>({sx, 0}), // A in cell
         }},
        {{
             Polynomial<1>({tx, 0}), // A in cell
             Polynomial<1>({0, 1}), // A <= B
             Polynomial<1>({sy, 0}), // A above valley
         }}
    }.translate_xy(-sx, -sx));

    // 3.
    // Via valley
    // A below valley
    // B above valley
    costs.push_back(ConstrainedBivariatePolynomial<2>{
        BivariatePolynomial<2>({{
            {{(sy * sy) / 2 + (ty * ty) / 2, -ty, 1./2}},
            {{-sy, 0, 0}},
            {{1./2, 0, 0}},
        }}),
        // B in cell + B above valley
        {std::clamp(sy, sx, tx), std::clamp(ty, sx, tx)},
        {{
             Polynomial<1>({sx, 0}), // A in cell
             Polynomial<1>({sy, 0}), // A below valley
         }},
        {{
             Polynomial<1>({tx, 0}), // A in cell
             Polynomial<1>({0, 1}), // A <= B
         }}
    }.translate_xy(-sx, -sx));

    // 4.
    // Via valley
    // A above valley
    // B below valley
    costs.push_back(ConstrainedBivariatePolynomial<2>{
        BivariatePolynomial<2>({{
            {{(sy * sy) / 2 + (ty * ty) / 2, -ty, 1./2}},
            {{-sy, 0, 0}},
            {{1./2, 0, 0}},
        }}),
        // B in cell + B below valley
        {std::clamp(ty, sx, tx), tx},
        {{
             Polynomial<1>({sx, 0}), // A in cell
         }},
        {{
             Polynomial<1>({tx, 0}), // A in cell
             Polynomial<1>({0, 1}), // A <= B
             Polynomial<1>({sy, 0}), // A above valley
         }}
    }.translate_xy(-sx, -sx));

    // 5.
    // Via valley
    // A below valley
    // B below valley
    costs.push_back(ConstrainedBivariatePolynomial<2>{
        BivariatePolynomial<2>({{
            {{(sy * sy) / 2 + (ty * ty) / 2, -ty, 1./2}},
            {{-sy, 0, 0}},
            {{1./2, 0, 0}},
        }}),
        // B in cell + B below valley
        {std::clamp(ty, sx, tx), tx},
        {{
             Polynomial<1>({sx, 0}), // A in cell
             Polynomial<1>({sy, 0}), // A below valley
         }},
        {{
             Polynomial<1>({tx, 0}), // A in cell
             Polynomial<1>({0, 1}), // A <= B
             Polynomial<1>({ty, 0}), // Via valley
         }}
    }.translate_xy(-sx, -sx));

    // 6.
    // Not via valley
    // A below valley
    // B below valley
    costs.push_back(ConstrainedBivariatePolynomial<2>{
        BivariatePolynomial<2>({{
            {{(sy * sy) / 2 - (ty * ty) / 2, -ty, 1./2}},
            {{-sy + 2 * ty, 0, 0}},
            {{-1./2, 0, 0}},
        }}),
        // B in cell + B below valley
        {std::clamp(ty, sx, tx), tx},
        {{
             Polynomial<1>({sx, 0}), // A in cell
             Polynomial<1>({ty, 0}), // A below valley + Not via valley
         }},
        {{
             Polynomial<1>({tx, 0}), // A in cell
             Polynomial<1>({0, 1}), // A <= B
         }}
    }.translate_xy(-sx, -sx));

    return costs;
}

#endif //TRAJECTORY_CLUSTERING_1D_L1_L1_H
