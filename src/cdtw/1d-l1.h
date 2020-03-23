#ifndef TRAJECTORY_CLUSTERING_1D_L1_H
#define TRAJECTORY_CLUSTERING_1D_L1_H

#include "cdtw.h"
#include "BivariatePolynomial.h"

PiecewisePolynomial<2> bottom_to_right_1(double sx, double sy, double tx, double ty) {
    // Not via valley
    // X below valley
    // Y below valley

    BivariatePolynomial<2> h({{
        {{(sy * sy) / 2 + (tx * tx) / 2, -tx, -1./2}},
        {{-sy, 2, 0}},
        {{-1./2, 0, 0}},
    }});

    // Y in cell + Y below valley
    Interval y_interval{sy, std::clamp(tx, sy, ty)};

    std::vector<Polynomial<1>> left_constraints;
    std::vector<Polynomial<1>> right_constraints;

    // X in cell
    left_constraints.push_back(Polynomial<1>({sx, 0}));
    right_constraints.push_back(Polynomial<1>({tx, 0}));

    // X below valley
    left_constraints.push_back(Polynomial<1>({sy, 0}));

    // Not via valley
    left_constraints.push_back(Polynomial<1>({0, 1}));

    return find_minimum(h, y_interval, left_constraints, right_constraints);
}

PiecewisePolynomial<2> bottom_to_right_2(double sx, double sy, double tx, double ty) {
    // Via valley
    // X below valley
    // Y below valley

    BivariatePolynomial<2> h({{
        {{(sy * sy) / 2 + (tx * tx) / 2, -tx, 1./2}},
        {{-sy, 0, 0}},
        {{1./2, 0, 0}},
    }});

    // Y in cell + Y below valley
    Interval y_interval{sy, std::clamp(tx, sy, ty)};

    std::vector<Polynomial<1>> left_constraints;
    std::vector<Polynomial<1>> right_constraints;

    // X in cell
    left_constraints.push_back(Polynomial<1>({sx, 0}));
    right_constraints.push_back(Polynomial<1>({tx, 0}));

    // X below valley
    left_constraints.push_back(Polynomial<1>({sy, 0}));

    // Via valley
    right_constraints.push_back(Polynomial<1>({0, 1}));

    return find_minimum(h, y_interval, left_constraints, right_constraints);
}

PiecewisePolynomial<2> bottom_to_right_3(double sx, double sy, double tx, double ty) {
    // Via valley
    // X above valley
    // Y below valley

    BivariatePolynomial<2> h({{
        {{(sy * sy) / 2 + (tx * tx) / 2, -tx, 1./2}},
        {{-sy, 0, 0}},
        {{1./2, 0, 0}},
    }});

    // Y in cell + Y below valley
    Interval y_interval{sy, std::clamp(tx, sy, ty)};

    std::vector<Polynomial<1>> left_constraints;
    std::vector<Polynomial<1>> right_constraints;

    // X in cell
    left_constraints.push_back(Polynomial<1>({sx, 0}));
    right_constraints.push_back(Polynomial<1>({tx, 0}));

    // X below valley
    left_constraints.push_back(Polynomial<1>({sy, 0}));

    return find_minimum(h, y_interval, left_constraints, right_constraints);
}

PiecewisePolynomial<2> bottom_to_right_4(double sx, double sy, double tx, double ty) {
    // Via valley
    // X below valley
    // Y above valley

    BivariatePolynomial<2> h({{
        {{(sy * sy) / 2 + (tx * tx) / 2, -tx, 1./2}},
        {{-sy, 0, 0}},
        {{1./2, 0, 0}},
    }});

    // Y in cell + Y above valley
    Interval y_interval{std::clamp(tx, sy, ty), ty};

    std::vector<Polynomial<1>> left_constraints;
    std::vector<Polynomial<1>> right_constraints;

    // X in cell
    left_constraints.push_back(Polynomial<1>({sx, 0}));
    right_constraints.push_back(Polynomial<1>({tx, 0}));

    // X below valley
    left_constraints.push_back(Polynomial<1>({sy, 0}));

    return find_minimum(h, y_interval, left_constraints, right_constraints);
}

PiecewisePolynomial<2> bottom_to_right_5(double sx, double sy, double tx, double ty) {
    // Via valley
    // X above valley
    // Y above valley

    // Via valley
    if (tx < sy) {
        return {};
    }

    BivariatePolynomial<2> h({{
        {{(sy * sy) / 2 + (tx * tx) / 2, -tx, 1./2}},
        {{-sy, 0, 0}},
        {{1./2, 0, 0}},
    }});

    // Y in cell + Y above valley
    Interval y_interval{std::clamp(tx, sy, ty), ty};

    std::vector<Polynomial<1>> left_constraints;
    std::vector<Polynomial<1>> right_constraints;

    // X in cell
    left_constraints.push_back(Polynomial<1>({sx, 0}));
    right_constraints.push_back(Polynomial<1>({tx, 0}));

    // X above valley
    right_constraints.push_back(Polynomial<1>({sy, 0}));

    return find_minimum(h, y_interval, left_constraints, right_constraints);
}

PiecewisePolynomial<2> bottom_to_right_6(double sx, double sy, double tx, double ty) {
    // Not via valley
    // X above valley
    // Y above valley

    // Not via valley
    if (tx > sy) {
        return {};
    }

    BivariatePolynomial<2> h({{
        {{2 * sy * tx - (sy * sy) / 2 - (tx * tx) / 2, -tx, 1./2}},
        {{-sy, 0, 0}},
        {{1./2, 0, 0}},
    }});

    // Y in cell + Y above valley
    Interval y_interval{std::clamp(tx, sy, ty), ty};

    std::vector<Polynomial<1>> left_constraints;
    std::vector<Polynomial<1>> right_constraints;

    // X in cell
    left_constraints.push_back(Polynomial<1>({sx, 0}));
    right_constraints.push_back(Polynomial<1>({tx, 0}));

    // X above valley
    right_constraints.push_back(Polynomial<1>({sy, 0}));

    return find_minimum(h, y_interval, left_constraints, right_constraints);
}

PiecewisePolynomial<2> bottom_to_right(double sx, double sy, double tx, double ty) {
    auto result = bottom_to_right_1(sx, sy, tx, ty);
//    std::cout << " ... " << result << '\n';
    fast_lower_envelope(result, bottom_to_right_2(sx, sy, tx, ty));
//    std::cout << " ... " << result << '\n';
    fast_lower_envelope(result, bottom_to_right_3(sx, sy, tx, ty));
//    std::cout << " ... " << result << '\n';
    fast_lower_envelope(result, bottom_to_right_4(sx, sy, tx, ty));
//    std::cout << " ... " << result << '\n';
    fast_lower_envelope(result, bottom_to_right_5(sx, sy, tx, ty));
//    std::cout << " ... " << result << '\n';
    fast_lower_envelope(result, bottom_to_right_6(sx, sy, tx, ty));
//    std::cout << " ... " << result << '\n';
    return result;
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
