#ifndef TRAJECTORY_CLUSTERING_2D_L2SQUARED_L1_H
#define TRAJECTORY_CLUSTERING_2D_L2SQUARED_L1_H

#include "cdtw.h"
#include "BivariatePolynomial.h"

//
// 2D + L2Squared image norm + L1 param norm
//

template<>
PiecewisePolynomial<3> CDTW<2, Norm::L2Squared, Norm::L1>::base_bottom(const Cell& cell) const {
    // Cell coordinates w.r.t. height function origin
    const double sx = cell.s.x - cell.mid.x;
    const double sy = cell.s.y - cell.mid.y;
    const double tx = cell.t.x - cell.mid.x;

    const double l = cell.l;

    return PiecewisePolynomial<3>({
        {{0, tx-sx}, Polynomial<3>({0, sx*sx + sy*sy + 2*l*sx*sy, sx + l*sy, 1./3})},
    });
}

template<>
std::vector<ConstrainedBivariatePolynomial<3>>
CDTW<2, Norm::L2Squared, Norm::L1>::bottom_to_right_costs(const Cell& cell) const {
    // Cell coordinates w.r.t. height function origin
    const double sx = cell.s.x - cell.mid.x;
    const double sy = cell.s.y - cell.mid.y;
    const double tx = cell.t.x - cell.mid.x;
    const double ty = cell.t.y - cell.mid.y;

    const double l = cell.l;

    std::vector<ConstrainedBivariatePolynomial<3>> costs;
    costs.reserve(5);

    // 1.
    // Not via valley
    // A below valley
    // B below valley
    costs.push_back(ConstrainedBivariatePolynomial<3>{
        BivariatePolynomial<3>({{
            {-((sx - tx)*(sx*sx + 3*l*sx*sy + 3*sy*sy + sx*tx + 3*l*sy*tx + tx*tx))/3.,-((-1 + l)*sx*sx) + 2*(-1 + l)*sx*sy + sy*sy + 2*sy*tx + l*tx*tx,(-1 + l)*sx + sy + tx,1./3},
            {-sx*sx - 2*l*sx*sy - sy*sy,-2*(-1 + l)*(sx - sy),-1 + l,0},
            {-sx - l*sy,1 - l,0,0},
            {-1./3,0,0,0},
        }}),
        // B in cell + B below valley
        Interval{0, std::clamp(tx-sy, 0., ty-sy)},
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
            {(-3*sx*sx*sy - sy*sy*sy + 3*sy*sy*tx + tx*tx*tx + l*(-sx*sx*sx - 3*sx*sy*sy + sy*sy*sy + 3*sy*tx*tx))/3.,2*sy*tx + l*(sy*sy + tx*tx),l*sy + tx,l/3.},
            {-2*sx*sy - l*(sx*sx + sy*sy),0,0,0},
            {-(l*sx) - sy,0,0,0},
            {-l/3.,0,0,0},
        }}),
        // B in cell + B below valley
        Interval{0, std::clamp(tx-sy, 0., ty-sy)},
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
            {-((sx - tx)*(sx*sx + 3*l*sx*sy + 3*sy*sy + sx*tx + 3*l*sy*tx + tx*tx))/3.,2*sy*tx + l*(sy*sy + tx*tx),l*sy + tx,l/3.},
            {-sx*sx - 2*l*sx*sy - sy*sy,0,0,0},
            {-sx - l*sy,0,0,0},
            {-1./3,0,0,0},
        }}),
        // B in cell + B below valley
        Interval{0, std::clamp(tx-sy, 0., ty-sy)},
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
            {-((sx - tx)*(3*sy*(sx + tx) + l*(sx*sx + 3*sy*sy + sx*tx + tx*tx)))/3.,sy*sy + 2*l*sy*tx + tx*tx,sy + l*tx,1./3},
            {-2*sx*sy - l*(sx*sx + sy*sy),0,0,0},
            {-(l*sx) - sy,0,0,0},
            {-l/3.,0,0,0},
        }}),
        // B in cell + B above valley
        Interval{std::clamp(tx-sy, 0., ty-sy), ty-sy},
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
                {(-sx*sx*sx - 3*l*sx*sx*sy - 3*sx*sy*sy - (-1 + l)*sy*sy*sy + 3*l*sy*sy*tx + 3*sy*tx*tx + l*tx*tx*tx)/3.,sy*sy + 2*l*sy*tx + tx*tx,sy + l*tx,1./3},
                {-sx*sx - 2*l*sx*sy - sy*sy,0,0,0},
                {-sx - l*sy,0,0,0},
                {-1./3,0,0,0},
            }}),
            // B in cell + B above valley
            Interval{std::clamp(tx-sy, 0., ty-sy), ty-sy},
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
                {-((sx - tx)*(sx*sx + 3*l*sx*sy + 3*sy*sy + sx*tx + 3*l*sy*tx + tx*tx))/3.,sy*sy + 2*l*sy*tx + tx*tx,sy + l*tx,1./3},
                {-sx*sx - 2*l*sx*sy - sy*sy,0,0,0},
                {-sx - l*sy,0,0,0},
                {-1./3,0,0,0},
            }}),
            // B in cell + B above valley
            Interval{std::clamp(tx-sy, 0., ty-sy), ty-sy},
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

template<>
std::vector<ConstrainedBivariatePolynomial<3>>
CDTW<2, Norm::L2Squared, Norm::L1>::bottom_to_top_costs(const Cell& cell) const {
    // Cell coordinates w.r.t. height function origin
    const double sx = cell.s.x - cell.mid.x;
    const double sy = cell.s.y - cell.mid.y;
    const double tx = cell.t.x - cell.mid.x;
    const double ty = cell.t.y - cell.mid.y;

    const double l = cell.l;

    std::vector<ConstrainedBivariatePolynomial<3>> costs;
    costs.reserve(6);

    // 1.
    // Not via valley
    // A above valley
    // B above valley
    costs.push_back(ConstrainedBivariatePolynomial<3>{
        BivariatePolynomial<3>({{
            {-((sy - ty)*(3*sx*sx + sy*sy + sy*ty + ty*ty + 3*l*sx*(sy + ty)))/3.,sx*sx - (-1 + l)*sy*sy + l*ty*ty + 2*sx*((-1 + l)*sy + ty),sx + (-1 + l)*sy + ty,1./3},
            {-sx*sx - 2*l*sx*sy - sy*sy,0,0,0},
            {-sx - l*sy,0,0,0},
            {-1./3,0,0,0},
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
            {((-1 + l)*sx*sx*sx - l*sy*sy*sy + ty*ty*ty + 3*sx*sx*(-(l*sy) + ty) - 3*sx*(sy*sy - l*ty*ty))/3.,2*sx*ty + l*(sx*sx + ty*ty),l*sx + ty,l/3.},
            {-sx*sx - 2*l*sx*sy - sy*sy,0,0,0},
            {-sx - l*sy,0,0,0},
            {-1./3,0,0,0},
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
            {-((sy - ty)*(3*sx*sx + sy*sy + sy*ty + ty*ty + 3*l*sx*(sy + ty)))/3.,2*sx*ty + l*(sx*sx + ty*ty),l*sx + ty,l/3.},
            {-2*sx*sy - l*(sx*sx + sy*sy),0,0,0},
            {-(l*sx) - sy,0,0,0},
            {-l/3.,0,0,0},
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
            {-((sy - ty)*(3*sx*(sy + ty) + l*(3*sx*sx + sy*sy + sy*ty + ty*ty)))/3.,sx*sx + 2*l*sx*ty + ty*ty,sx + l*ty,1./3},
            {-sx*sx - 2*l*sx*sy - sy*sy,0,0,0},
            {-sx - l*sy,0,0,0},
            {-1./3,0,0,0},
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
            {(-((-1 + l)*sx*sx*sx) - sy*sy*sy + l*ty*ty*ty - 3*sx*sx*(sy - l*ty) + 3*sx*(-(l*sy*sy) + ty*ty))/3.,sx*sx + 2*l*sx*ty + ty*ty,sx + l*ty,1./3},
            {-2*sx*sy - l*(sx*sx + sy*sy),0,0,0},
            {-(l*sx) - sy,0,0,0},
            {-l/3.,0,0,0},
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
            {-((sy - ty)*(3*sx*sx + sy*sy + sy*ty + ty*ty + 3*l*sx*(sy + ty)))/3.,sx*sx + 2*l*sx*ty + ty*ty,sx + l*ty,1./3},
            {-sx*sx - l*sy*sy + (-1 + l)*ty*ty - 2*sx*(sy + (-1 + l)*ty),0,0,0},
            {-sx - sy + ty - l*ty,0,0,0},
            {-1./3,0,0,0},
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

#endif //TRAJECTORY_CLUSTERING_2D_L2SQUARED_L1_H
