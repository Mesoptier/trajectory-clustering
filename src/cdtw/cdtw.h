#ifndef CDTW_H
#define CDTW_H
#define FAST_LOWER_ENVELOPE true

#include <vector>
#include <array>
#include <ostream>
#include <cmath>
#include <iostream>
#include <set>
#include <algorithm>
#include <iomanip>

#include "Interval.h"
#include "Polynomial.h"
#include "PiecewisePolynomial.h"
#include "BivariatePolynomial.h"
#include "ConstrainedBivariatePolynomial.h"
#include "naive_lower_envelope.h"
#include "fast_lower_envelope.h"
#include "../IntegralFrechet/Cell.h"


template<size_t D>
void verify_minimum(const PiecewisePolynomial<D>& f, const std::vector<ConstrainedBivariatePolynomial<D>>& hs);

/**
 * Given a bivariate polynomial function H(x,y) computes a univariate piecewise polynomial that returns for each Y
 * the minimum value bounded by the interval and left/right constraints (i.e. y maps to min_x H(x,y))
 *
 * @tparam D - Degree of the polynomials
 * @param h - Bivariate polynomial function H(x,y)
 * @param interval_y - Domain of Y
 * @param left_constraints - Constraints of shape a*y+b >= 0
 * @param right_constraints - Constraints of shape a*y+b <= 0
 * @return
 */
template<size_t D>
PiecewisePolynomial<D> find_minimum(
    const BivariatePolynomial<D>& h,
    const Interval_c& interval_y,
    std::vector<Polynomial<1>> left_constraints,
    std::vector<Polynomial<1>> right_constraints,
    const ConstrainedBivariatePolynomial<D>& f
);

template<size_t D>
void fast_lower_envelope(PiecewisePolynomial<D>& left, const PiecewisePolynomial<D>& right);

struct Entry
    {
        // Functions along x- and y-axis respectively
        PiecewisePolynomial<2> bottom;
        PiecewisePolynomial<2> left;
    };

template<size_t dimension, Norm image_norm, Norm param_norm>
class CDTW
{
public:
    static const size_t D =
        (image_norm == Norm::L1 && param_norm == Norm::L1) ? 2 :
        (image_norm == Norm::L2Squared && param_norm == Norm::L1) ? 3 :
        0;

    std::vector<Point> warping_path;

    void write_heat_map(const Curve& curve1, const Curve& curve2, std::string norm="L1");

    CDTW(Curve const& _curve1, Curve const& _curve2);

    double cost() const;

    void print_complexity() const;

    void output_visualization_data() const;

    const std::vector<std::vector<Entry>>& get_functions() const;

    void check_table();

    void find_max_pieces();

    int max_pieces;

    void write_warping_path();

private:
    const Curve& curve1;
    const Curve& curve2;

    size_t n;
    size_t m;

    double cdtw_cost;

    // Dynamic program
    std::vector<std::vector<Entry>> in_functions;

    // Result
    PiecewisePolynomial<D> out_top;
    PiecewisePolynomial<D> out_right;

    // Internal methods
    template<class Iterator>
    PiecewisePolynomial<2> propagate(
        const std::vector<ConstrainedBivariatePolynomial<D>>& cell_costs,
        Iterator pieces_it,
        Iterator pieces_end
    ) const;


    PiecewisePolynomial<2> bottom_to_right(const PiecewisePolynomial<D>& in, const Cell& cell) const;

    PiecewisePolynomial<2> bottom_to_top(const PiecewisePolynomial<D>& in, const Cell& cell) const;

    void get_path_segs(PathType path_type, Boundaries boundaries, double in, double out,
    std::vector<std::vector<double>>& output, int i, int j, Line axis, bool left);

    void compute_warping_path(std::vector<std::vector<double>>& output, int i, int j, double y, bool left);

    void process_cell(int i, int j, const Curve& curve1, const Curve& curve2);

    // Internal methods specific to each case
    PiecewisePolynomial<2> base_bottom(const Cell& cell) const;
    std::vector<ConstrainedBivariatePolynomial<2>> bottom_to_right_costs(const Cell& cell) const;
    std::vector<ConstrainedBivariatePolynomial<2>> bottom_to_top_costs(const Cell& cell) const;
    
};

// TODO: Remove temporary checks to find problematic cases for Sampson's proof.
template<size_t D>
void check_function(const PiecewisePolynomial<D>& f);

#endif

