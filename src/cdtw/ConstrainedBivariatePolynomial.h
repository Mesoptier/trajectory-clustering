#ifndef TRAJECTORY_CLUSTERING_CONSTRAINEDBIVARIATEPOLYNOMIAL_H
#define TRAJECTORY_CLUSTERING_CONSTRAINEDBIVARIATEPOLYNOMIAL_H

#include <ostream>
#include "PiecewisePolynomial.h"
#include "BivariatePolynomial.h"

// TODO: Rename to BivariatePolynomialPiece

template<size_t D>
struct ConstrainedBivariatePolynomial
{
    // f(x, y)
    BivariatePolynomial<D> f;
    // y is valid in this interval
    Interval y_interval;
    // x is valid when g(y) <= x <= h(y) for each g in left_constraints and h in right_constraints
    std::vector<Polynomial<1>> left_constraints;
    std::vector<Polynomial<1>> right_constraints;

    PathType path_type = NONE;
    Boundaries boundaries;
    Line axis;
    std::string history;
    double test_value = 0;
    double xmin = 0;
    double ymax = 0;

    inline static int count = 0;
    int _id;

    void transpose() {

        auto dir = axis.direction;
        auto p = Point(0, axis.getY(0));

        axis = Line::fromTwoPoints(
            Point(p.y, p.x), Point(p.y + dir.y, p.x + dir.x)
        );

        switch (boundaries) {
            case BR: {
                boundaries = LT;
                break;
            }
            case BT: {
                boundaries = LR;
                break;
            }
            case LR: {
                boundaries = BT;
                break;
            }
            case LT: {
                boundaries = BR;
                break;
            }
        }

        switch (path_type) {
            case UA: {
                path_type = AU;
                break;
            }
            case AU: {
                path_type = UA;
                break;
            }
            case UAA: {
                path_type = AAU;
                break;
            }
            case AAU: {
                path_type = UAA;
                break;
            }
            case AAA: {
                path_type = UAU;
                break;
            }
            case UAU: {
                path_type = AAA;
                break;
            }
        }
    }
    
    /**
     * Sum this function with a function in x.
     * Returns f(x, y) + g(x) with updated constraints.
     */
    ConstrainedBivariatePolynomial<D> add_x(const PolynomialPiece<D>& g) const;

    ConstrainedBivariatePolynomial<D> translate_xy(double cx, double cy) const {
        auto result = *this;
        result.f = result.f.translate_xy(cx, cy);
        result.y_interval.min += cy;
        result.y_interval.max += cy;
        for (auto& c : result.left_constraints) {
            c = c.translate_xy(cy, cx);
        }
        for (auto& c : result.right_constraints) {
            c = c.translate_xy(cy, cx);
        }
        return result;
    }

    ConstrainedBivariatePolynomial<D> multiply(double a) const {
        auto result = *this;
        result.f *= a;
        return result;
    }

    PolynomialPiece<D> slice_at_y(double c) const {
        return {interval_at_y(c), f.slice_at_y(c)};
    }

    Interval interval_at_y(double y) const {
        double min = -std::numeric_limits<double>::infinity();
        double max = std::numeric_limits<double>::infinity();

        for (const auto& lc : left_constraints) {
            min = std::max(min, lc(y));
        }
        for (const auto& rc : right_constraints) {
            max = std::min(max, rc(y));
        }

        return { min, max };
    }

    bool constraints_feasible() {
        
    }

    friend std::ostream& operator<<(std::ostream& os, const ConstrainedBivariatePolynomial& polynomial) {
        os << std::fixed;
        os << "{(" << polynomial.f << "), ";
        os << "((" << polynomial.y_interval << ") /. x->y)";
        for (auto& c : polynomial.left_constraints) {
            os << " && ((" << c << ") /. x->y) <= x";
        }
        for (auto& c : polynomial.right_constraints) {
            os << " && ((" << c << ") /. x->y) >= x";
        }
        os << "}";
        return os;
    }

    double operator()(Point p) {
        return f(p);
    }

    std::vector<Polynomial<1>> remove_redundant_const(
        std::vector<Polynomial<1>> init_set, std::string side) {

        std::vector<Polynomial<1>> result;

        if (approx_equal(init_set[0].coefficients[0], 1.4142135623730951) && approx_equal(init_set[0].coefficients[1], -1.0000049999625003))
            std::cout << "...\n";

        for (auto c: init_set) {  

            bool undominated = true;
            double c_max = std::max(c(y_interval.min), c(y_interval.max));
            double c_min = std::min(c(y_interval.min), c(y_interval.max));
            
            for (auto oth_c: init_set) {

                double oth_c_max = std::max(oth_c(y_interval.min), oth_c(y_interval.max));
                double oth_c_min = std::min(oth_c(y_interval.min), oth_c(y_interval.max));

                if (c != oth_c)
                    undominated &= side == "right" ? (c_min < oth_c_min ||  c_max < oth_c_max)
                    : (c_min > oth_c_min ||  c_max > oth_c_max);
            }

            if (undominated) {
                if (std::find(result.begin(), result.end(), c) == result.end())
                    result.push_back(c);
            }
        }

        if (result.size() == 0)
            std::cout << "...\n";

        return result;
    }

    void clean_constraints() {
        // auto new_left = std::vector<Polynomial<1>>();

        // bool zero_seen = false;

        // for (auto c: left_constraints) 
        //     if (c.coefficients[0] >= 0 || c.coefficients[1] >= 0) {
        //         if (c.coefficients[0] == 0 && c.coefficients[1] == 0 && !zero_seen) {
        //             new_left.push_back(c);
        //             zero_seen = true;
        //         } else if (c.coefficients[0] > 0 || c.coefficients[1] > 0)
        //             new_left.push_back(c);
        //     }

  

        // left_constraints = new_left;

        left_constraints = remove_redundant_const(left_constraints, "left");
        right_constraints = remove_redundant_const(right_constraints, "right");
    }

    friend ConstrainedBivariatePolynomial operator+(
        ConstrainedBivariatePolynomial& l, ConstrainedBivariatePolynomial& r
    ) {
        auto polynomial = l.f + r.f;
        auto y_interval = l.y_interval.intersect(r.y_interval);

        auto left_constraints = std::vector<Polynomial<1>>();
        auto right_constraints = std::vector<Polynomial<1>>();

        for (auto& poly: l.left_constraints) 
            left_constraints.push_back(poly);

        for (auto& poly: r.left_constraints) {
            bool new_constraint = true;
            for (auto& other: left_constraints)
                if (other == poly)
                    new_constraint = false;
            if (new_constraint)
                left_constraints.push_back(poly);
        }

        for (auto& poly: l.right_constraints)
            right_constraints.push_back(poly);

        for (auto& poly: r.right_constraints) {
            bool new_constraint = true;
            for (auto& other: right_constraints)
                if (other == poly)
                    new_constraint = false;
            if (new_constraint)
                right_constraints.push_back(poly);
        }

        auto output = ConstrainedBivariatePolynomial{
            polynomial,
            y_interval,
            left_constraints,
            right_constraints
        };

        output.history = l.history + "_" + r.history;
    

        return output;
    }
};

template<size_t D>
ConstrainedBivariatePolynomial<D> ConstrainedBivariatePolynomial<D>::add_x(const PolynomialPiece<D>& g) const {
    auto result = *this;
    result.f = result.f.add_x(g.polynomial);
    result.left_constraints.push_back(Polynomial<1>({g.interval.min, 0}));
    result.right_constraints.push_back(Polynomial<1>({g.interval.max, 0}));
    return result;
}

#endif //TRAJECTORY_CLUSTERING_CONSTRAINEDBIVARIATEPOLYNOMIAL_H
