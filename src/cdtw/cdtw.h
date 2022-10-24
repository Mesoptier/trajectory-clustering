#ifndef CDTW_H
#define CDTW_H

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

    void write_heat_map(Curve& new_curve1, Curve& new_curve2, std::string norm="L1");

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
    std::vector<std::vector<double>>& output, Curve& c1, Curve& c2, int i, int j, Line axis, bool left);

    void compute_warping_path(Curve& new_curve1, Curve& new_curve2, 
    std::vector<std::vector<double>>& output, int i, int j, double y, bool left);

    void process_cell(int i, int j, Curve& curve1, Curve& curve2);

    // Internal methods specific to each case
    PiecewisePolynomial<2> base_bottom(const Cell& cell) const;
    std::vector<ConstrainedBivariatePolynomial<2>> bottom_to_right_costs(const Cell& cell) const;
    std::vector<ConstrainedBivariatePolynomial<2>> bottom_to_top_costs(const Cell& cell) const;

public:
    CDTW(Curve const& _curve1, Curve const& _curve2);

    double cost() const;

    void print_complexity() const;

    void output_visualization_data() const;

    const std::vector<std::vector<Entry>>& get_functions() const;

    void check_table();

    void find_max_pieces();

    int max_pieces;
};

// TODO: Remove temporary checks to find problematic cases for Sampson's proof.
template<size_t D>
void check_function(const PiecewisePolynomial<D>& f);

#endif
// template<size_t dimension, Norm image_norm, Norm param_norm>
// CDTW<dimension, image_norm, param_norm>::CDTW(Curve& _curve1, Curve& _curve2) :
//     curve1(_curve1), curve2(_curve2),
//     n(curve1.size()), m(curve2.size()),
//     in_functions(n, std::vector<Entry>(m))
// {
//     struct Function {
//         PiecewisePolynomial<D> f;
//         distance_t min_cost;
//         bool open; // Whether this function is currently in the open_set
//     };

//     // Curve new_curve1("", curve1.get_points());
//     // Curve new_curve2("", curve2.get_points());

//     // if (dimension == 2 && image_norm == Norm::L1 && param_norm == Norm::L1) {

//         int new_n;
//         int new_m;

//         std::vector<CPoint> new_c_points_1 = std::vector<CPoint>();
//         std::vector<CPoint> new_c_points_2 = std::vector<CPoint>();

//         for (int i = 0; i < curve1.size(); ++i) {
//             new_c_points_1.push_back(CPoint(i, 0));
//         }

//         for (int i = 0; i < curve2.size(); ++i) {
//             new_c_points_2.push_back(CPoint(i, 0));
//         }


//         for (int i = 0; i < curve1.size()-1; ++i) {
//             auto s1 = curve1[i];
//             auto t1 = curve1[i+1];

//             const auto l1 = Line::fromTwoPoints(s1, t1);

//             for (int j = 0; j < curve2.size()-1; ++j) {

//                 auto s2 = curve2[j];
//                 auto t2 = curve2[j+1];

//                 const auto l2 = Line::fromTwoPoints(s2, t2);

//                 if (!isParallel(l1, l2)) {
//                     const auto p = intersect(l1, l2);
//                     Point mid = p;

//                     distance_t c1_dist = s1.dist(mid);
//                     distance_t c2_dist = s2.dist(mid);

//                     if ((mid.x - s1.x)*(mid.x - t1.x) < 0
//                     && (mid.x - s2.x)*(mid.x - t2.x) < 0) {
//                         std::cout << c1_dist/(s1.dist(t1)) << std::endl;
//                         std::cout << c2_dist/(s2.dist(t2)) << std::endl;
//                         // std::cout << s1 << std::endl;
//                         // std::cout << t1 << std::endl;
//                         // std::cout << s2 << std::endl;
//                         // std::cout << t2 << std::endl;
//                         // std::cout << mid << std::endl;
//                         new_c_points_1.push_back(CPoint(i, c1_dist/(s1.dist(t1))));
//                         new_c_points_2.push_back(CPoint(j, c2_dist/(s2.dist(t2))));
//                     }
//                 }
//             }
//         }

//         std::sort(new_c_points_1.begin(), new_c_points_1.end());
//         std::sort(new_c_points_2.begin(), new_c_points_2.end());
    
//         Points new_points_1 = Points();

//         for (auto c: new_c_points_1)
//             new_points_1.push_back(curve1.interpolate_at(c));

//         Points new_points_2 = Points();

//         for (auto c: new_c_points_2)
//             new_points_2.push_back(curve2.interpolate_at(c));

//         Curve new_curve1 = Curve("", new_points_1);
//         Curve new_curve2 = Curve("", new_points_2);
//     // }

//     std::unordered_map<Coord, Function, typename Coord::hash> functions;
//     std::priority_queue<PQNode, std::vector<PQNode>, std::greater<>> open_set;

//     {
//         const Cell cell(new_curve1[0], new_curve2[0], new_curve1[1], new_curve2[1]);

//         const Function f1{base_bottom(cell), 0, true};
//         const Coord c1{0, 0, Side::BOTTOM};
//         functions.insert_or_assign(c1, f1);
//         open_set.push({f1.f.min_value(), c1});
//     }

//     distance_t min_result_cost = std::numeric_limits<distance_t>::infinity();

//     size_t nodes_handled = 0;
//     size_t nodes_skipped = 0;
//     size_t nodes_opened = open_set.size();
//     size_t nodes_reopened = 0;

//     while (!open_set.empty()) {
//         const PQNode node = open_set.top();
//         open_set.pop();

//         std::cout << node.coord << std::endl;

// //        std::cout << min_result_cost << " " << node << '\n';
//         const Coord coord = node.coord;
//         Function& f = functions.at(coord);
//         f.open = false;

//         // We can no longer find a lower min_result_cost, so the search is finished
//         if (node.min_cost >= min_result_cost) {
//             break;
//         }

//         if (node.min_cost > f.min_cost) {
//             ++nodes_skipped;
//             continue;
//         }

//         if (coord.i1 == new_curve1.size() - 1 || coord.i2 == new_curve2.size() - 1) {
//             if (coord.i1 == new_curve1.size() - 2 || coord.i2 == new_curve2.size() - 2) {
//                 min_result_cost = std::min(min_result_cost, f.f.right_value());
//             }
//             ++nodes_skipped;
//             continue;
//         }

//         ++nodes_handled;

//         // Cell and transposed cell
//         const Cell cell(new_curve1[coord.i1], new_curve2[coord.i2], new_curve1[coord.i1 + 1], new_curve2[coord.i2 + 1]);
//         const Cell cell_t(new_curve2[coord.i2], new_curve1[coord.i1], new_curve2[coord.i2 + 1], new_curve1[coord.i1 + 1]);

//         for (size_t i = 0; i < 2; ++i) { // 0 -> right, 1 -> top
//             const PiecewisePolynomial<D> out_f = i == 0
//                 ? (coord.side == Side::BOTTOM ? bottom_to_right(f.f, cell) : bottom_to_top(f.f, cell_t))
//                 : (coord.side == Side::BOTTOM ? bottom_to_top(f.f, cell) : bottom_to_right(f.f, cell_t));
//             const Coord out_coord = i == 0
//                 ? coord.right()
//                 : coord.top();

            

//             distance_t min_cost = out_f.min_value();

//             auto it = functions.find(out_coord);
//             if (it == functions.end()) {
//                 functions.emplace(out_coord, Function{out_f, min_cost, true});
//                 open_set.push({min_cost, out_coord});
//                 ++nodes_opened;
//             } else {
//                 Function& prev_out = it->second;
//                 // TODO: Only re-open function if lower_envelope has changed
//                 fast_lower_envelope(prev_out.f, out_f);

//                 if (min_cost < prev_out.min_cost) {
//                     prev_out.open = true;
//                     open_set.push({min_cost, out_coord});
//                     ++nodes_reopened;
//                 } else { // min_cost == prev_out.min_cost
//                     if (!prev_out.open) {
//                         prev_out.open = true;
//                         open_set.push({min_cost, out_coord});
//                         ++nodes_reopened;
//                     }
//                 }
//             }
//         }
//     }


//     cdtw_cost = min_result_cost;
//     std::cout << min_result_cost << "\n";
//     std::cout << "nodes_handled: " << nodes_handled << "\n";
//     std::cout << "nodes_skipped: " << nodes_skipped << "\n";
//     std::cout << "nodes_opened: " << nodes_opened << "\n";
//     std::cout << "nodes_reopened: " << nodes_reopened << "\n";
//     std::cout << "nodes remaining: " << open_set.size() << "\n";
// }
