#include <iostream>
#include <chrono>
#include <fstream>
#include <random>

#include "src/io.h"
#include "src/Curve.h"
#include "src/IntegralFrechet/IntegralFrechet.h"
#include "src/SymmetricMatrix.h"
#include "src/clustering/pam.h"
#include "src/IntegralFrechet/MatchingBand.h"
#include "src/simplification/imaiiri.h"
#include "src/DTW/dtw.h"
#include "src/cdtw/cdtw.h"
#include "src/cdtw/1d-l1-l1.h"
#include "src/cdtw/2d-l1-l1.h"
#include "src/cdtw/1d-l2squared-l1.h"
#include "src/cdtw/2d-l2squared-l1.h"

// namespace {
//
// I/O helpers
//

// std::vector<Curve> read_curves(const std::string& filename) {
//     auto start_time = std::chrono::high_resolution_clock::now();

//     const auto curves = io::read_curves(filename);

//     auto end_time = std::chrono::high_resolution_clock::now();
//     auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
//     std::cout << "read_curves: " << duration << "ms\n";

//     return curves;
// }

// SymmetricMatrix read_matrix(const std::string& filename) {
//     auto start_time = std::chrono::high_resolution_clock::now();

//     std::ifstream file(filename);
//     const auto matrix = SymmetricMatrix::read(file);

//     auto end_time = std::chrono::high_resolution_clock::now();
//     auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
//     std::cout << "read_matrix: " << duration << "ms\n";

//     return matrix;
// }

// void export_matrix(const SymmetricMatrix& matrix, const std::string& filename) {
//     std::ofstream file(filename);
//     matrix.write(file);
//     file.close();
// }


// //
// // Algorithms
// //

// SymmetricMatrix compute_distance_matrix(const std::vector<Curve>& curves) {
//     size_t n = curves.size();
//     SymmetricMatrix distance_matrix(n);

//     auto start_time = std::chrono::high_resolution_clock::now();

//     for (size_t i = 0; i < n; ++i) {
//         const auto curve1 = curves[i];

//         // Don't compute the j == i case, since the cost is always 0 (comparing identical curves)
//         distance_matrix.at(i, i) = 0;

//         for (size_t j = i + 1; j < n; ++j) {
//             const auto curve2 = curves[j];

//             IntegralFrechet alg(curve1, curve2, ParamMetric::LInfinity_NoShortcuts, 100, nullptr);
//             const auto result = alg.compute_matching();

//             // TODO: Weigh by max or sum?
//             // TODO: Move to IntegralFrechet class
//             distance_t averaged_cost = result.cost / std::max(curve1.curve_length(), curve2.curve_length());
//             distance_matrix.at(i, j) = averaged_cost;
//         }

//         std::cout << i << '\n';
//     }

//     auto end_time = std::chrono::high_resolution_clock::now();
//     auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
//     std::cout << "compute_distance_matrix: " << duration << "ms\n";

//     return distance_matrix;
// }

// void compute_clusters(const SymmetricMatrix& distance_matrix, size_t k) {
//     auto start_time = std::chrono::high_resolution_clock::now();

//     clustering::pam::compute(distance_matrix.n, k, distance_matrix);

//     auto end_time = std::chrono::high_resolution_clock::now();
//     auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
//     std::cout << "compute_clusters: " << duration << "ms\n";
// }

// IntegralFrechet::MatchingResult compute_matching(const Curve& curve1, const Curve& curve2, distance_t resolution, const MatchingBand* const band) {
//     auto start_time = std::chrono::high_resolution_clock::now();

//     IntegralFrechet alg(curve1, curve2, ParamMetric::LInfinity_NoShortcuts, resolution, band);
//     const auto result = alg.compute_matching();

//     auto end_time = std::chrono::high_resolution_clock::now();
//     auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
//     std::cout << "compute_matching: " << duration << "ms\n";

//     return result;
// }

// MatchingBand compute_band(const Curve& curve1, const Curve& curve2, const Points& matching, distance_t radius) {
//     auto start_time = std::chrono::high_resolution_clock::now();

//     const MatchingBand band(curve1, curve2, matching, radius);

//     auto end_time = std::chrono::high_resolution_clock::now();
//     auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
//     std::cout << "compute_band: " << duration << "ms\n";

//     return band;
// }

// //
// // Experiments
// //

// void experiment_with_or_without_bands() {
//     struct Stat {
//         distance_t band_cost;
//         long long int band_time;
//         distance_t no_band_cost;
//         long long int no_band_time;
//     };
//     std::vector<Stat> stats;

//     bool maintain_lengths = true;
//     distance_t simple_resolution = 10;
//     distance_t resolution = 1;
//     distance_t band_radius = 1;

//     const auto curves = read_curves("data/characters/data");

//     std::random_device rd;
//     std::mt19937_64 generator(rd());
//     std::uniform_int_distribution<size_t> distribution(0, curves.size() - 1);

//     size_t k_max = 1000;
//     for (size_t k = 0; k < k_max; ++k) {
//         auto i = distribution(generator);
//         auto j = distribution(generator);

//         const auto curve1 = curves.at(i);
//         const auto curve2 = curves.at(j);

//         // With band
//         auto start_band = std::chrono::high_resolution_clock::now();
//         const auto result_alt = compute_matching(curve1.simplify(maintain_lengths), curve2.simplify(maintain_lengths), simple_resolution,nullptr);
//         const auto band = compute_band(curve1, curve2, result_alt.matching, band_radius);
//         const auto result_band = compute_matching(curve1, curve2, resolution, &band);
//         auto end_band = std::chrono::high_resolution_clock::now();
//         auto duration_band = std::chrono::duration_cast<std::chrono::milliseconds>(end_band - start_band).count();

//         // Without band
//         auto start_no_band = std::chrono::high_resolution_clock::now();
//         const auto result_no_band = compute_matching(curve1, curve2, resolution, nullptr);
//         auto end_no_band = std::chrono::high_resolution_clock::now();
//         auto duration_no_band = std::chrono::duration_cast<std::chrono::milliseconds>(end_no_band - start_no_band).count();

//         stats.push_back({
//             result_band.cost,
//             duration_band,
//             result_no_band.cost,
//             duration_no_band,
//         });

//         std::cout << (k + 1) << " / " << k_max << '\n';
//     }

//     std::ofstream file("data/out/band_stats.tsv");
//     file << "band_cost\tband_time\tno_band_cost\tno_band_time\n";
//     for (const auto& stat : stats) {
//         file << stat.band_cost << '\t' << stat.band_time << '\t' << stat.no_band_cost << '\t' << stat.no_band_time << '\n';
//     }
//     file.close();
// }

// void experiment_visualize_band() {
//     bool maintain_lengths = true;

//     const auto curve1 = io::read_curve("data/characters/data/a0001.txt");
//     const auto curve2 = io::read_curve("data/characters/data/a0002.txt");

//     const auto result_no_band = compute_matching(curve1, curve2, 1, nullptr);
//     io::export_points("data/out/debug_points_old.csv", result_no_band.search_stat.nodes_as_points);
//     io::export_points("data/out/matching3.csv", result_no_band.matching);
//     std::cout << "matching (no band) cost: " << result_no_band.cost << '\n';

//     auto start_time = std::chrono::high_resolution_clock::now();

//     const auto result_alt = compute_matching(curve1.simplify(maintain_lengths), curve2.simplify(maintain_lengths), 10,nullptr);
//     const auto band = compute_band(curve1, curve2, result_alt.matching,     2);
//     const auto result = compute_matching(curve1, curve2, 1, &band);

//     auto end_time = std::chrono::high_resolution_clock::now();
//     auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
//     std::cout << "main: " << duration << "ms\n";

//     io::export_points("data/out/curve1.csv", curve1.get_points());
//     io::export_points("data/out/curve2.csv", curve2.get_points());
//     io::export_points("data/out/matching1.csv", result.matching);
//     io::export_points("data/out/matching2.csv", result_alt.matching);
//     io::export_points("data/out/debug_points.csv", result.search_stat.nodes_as_points);

//     std::cout << "matching cost: " << result.cost << '\n';
//     std::cout << "matching (alt) cost: " << result_alt.cost << '\n';
// }

// void experiment_compare_heuristic_vs_extact_cdtw() {
// //    const std::string data_dir = "D:/TUe/MSc Project/code/data";
//     const std::string data_dir = "data";

//     auto curve1 = io::read_curve(data_dir + "/characters/data/a0001.txt");
//     auto curve2 = io::read_curve(data_dir + "/characters/data/a0002.txt");

// //    const auto curve1 = io::read_curve("data/characters/data/a0003.txt").simplify(false);
// //    const auto curve2 = io::read_curve("data/characters/data/a0002.txt").simplify(false);

// //    const Curve curve2("segment", {{0, 0}, {10, 10}});

// //    io::export_points("data/out/curve1.csv", curve1.get_points());
// //    io::export_points("data/out/curve2.csv", curve2.get_points());


//     // Heuristic
//     {
//         auto start_time = std::chrono::high_resolution_clock::now();

//         IntegralFrechet heuristic_alg(curve1, curve2, ParamMetric::L1, .1);
//         const auto heuristic_res = heuristic_alg.compute_matching();

//         auto end_time = std::chrono::high_resolution_clock::now();
//         auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
//         std::cout << "heuristic cost: " << heuristic_res.cost << '\n';
//         std::cout << "heuristic time: " << duration << "ms\n";
//         std::cout << "heuristic nodes_handled: " << heuristic_res.search_stat.nodes_handled << "\n";
//         std::cout << "heuristic nodes_skipped: " << heuristic_res.search_stat.nodes_skipped << "\n";
//         std::cout << "heuristic nodes_opened: " << heuristic_res.search_stat.nodes_opened << "\n";
//     }

//     // Exact
//     {
//         auto start_time = std::chrono::high_resolution_clock::now();

//         CDTW<2, Norm::L2Squared, Norm::L1> exact_alg(curve1, curve2);

//         auto end_time = std::chrono::high_resolution_clock::now();
//         auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
// //        std::cout << "exact cost: " << exact_alg.cost() << '\n';
//         std::cout << "exact time: " << duration << "ms\n";

// //        exact_alg.output_visualization_data();
// //        exact_alg.print_complexity();
//     }

// }
// }

// void exact_2dl1l1_heuristic_comp() {

//     using _CDTW = CDTW<2, Norm::L1, Norm::L1>;

//     auto characters = io::read_curves("data/characters/data");
//     auto test_curves = std::vector<Curve>();

//     for (int i = 0; i < 20; ++i) {
//         test_curves.push_back(characters[10*i].slice(1, 20));
//     }


//     double max = -1;
//     double min = 10000000;
//     double sum = 0;
    
//     for (int i = 0; i < test_curves.size(); ++i) 
//         for (int j = i+1; j < test_curves.size(); ++j) {
//             if (j == 6 && i == 0)
//                 std::cout << "ahh\n";
                
//             auto cdtw = _CDTW(test_curves[i], test_curves[j]);
//             IntegralFrechet heuristic_alg(test_curves[i], test_curves[j], ParamMetric::L1, .1);
//             const auto heuristic_res = heuristic_alg.compute_matching();

//             std::cout << "exact: " << cdtw.cost() << std::endl;
//             std::cout << "heur: " << heuristic_res.cost << std::endl;
//             std::cout << i << std::endl;
//             std::cout << j << std::endl;

//             double ratio = cdtw.cost() / heuristic_res.cost;
//             std::cout << ratio << std::endl;
//             max = std::max(max, ratio);
//             min = std::min(min, ratio);
//             sum += ratio;

//         };

//         double average_ratio = sum / (test_curves.size()*(test_curves.size()-1) / 2);

//         std::cout << "average ratio: " << average_ratio << std::endl;
//         std::cout << "max ratio: " << max << std::endl;
//         std::cout << "min ratio: " << min << std::endl;
    
// }

void test_case_1() {
    Point p1(0, 0);
    Point p2(0, 1);
    Point p3(0, 2);
    Point p4(0, 3);

    Point q1(0 ,0);
    Point q2(1, 1);
    Point q3(0, 2);
    Point q4(1, 3);

    auto c1 = Curve("", {
        p1, p2, p3, p4
    });

    auto c2 = Curve("", {
        q1, q2, q3, q4
    });

    using _CDTW = CDTW<2, Norm::L1, Norm::L1>;
    auto cdtw = _CDTW(c1, c2);

    std::cout << "test case 1 cost: " << cdtw.cost() << std::endl;
}

void test_case_2() {

    Point p1(0, 0);
    Point p2(1, 1);
    Point p3(6, 2);
    Point p4(7, 3);
    Point p5(1, 0);
    Point p6(0, 1);

    Point q1(0, 0);
    Point q2(1, 1);
    Point q3(6, 2);
    Point q4(7, 3);
    Point q5(1, 0);
    Point q6(0, 1);

    Curve c1("", {
        p1, p2, p3, p4, p5, p6
    });

    Curve c2("", {
        q1, q2, q3, q4, q5, q6
    });

    using _CDTW = CDTW<2, Norm::L1, Norm::L1>;
    auto cdtw = _CDTW(c1, c2);

    std::cout << "test case 2 cost: " << cdtw.cost() << std::endl;
}

void test_case_3() {
    Point p1(0, 0);
    Point p2(0, 0.25);
    Point p3(0, 0.5);
    Point p4(0, 0.75);
    Point p5(0, 1);

    Point q1(0, 0);
    Point q2(0.25, 0.25);
    Point q3(0.5, 0.5);
    Point q4(0.75, 0.75);
    Point q5(1, 1);

    Curve c1("", {
        p1, p2, p3, p4, p5
    });

    Curve c2("", {
        q1, q2, q3, q4, q5
    });

    using _CDTW = CDTW<2, Norm::L1, Norm::L1>;
    auto cdtw = _CDTW(c1, c2);

    std::cout << "test case 3 cost: " << cdtw.cost() << std::endl;
}

void test_case_4() {
    Point p1(0, 0);
    Point p2(0, 0.2);
    Point p3(0, 0.4);
    Point p4(0, 0.6);
    Point p5(0, 0.8);
    Point p6(0, 1);

    Point q1(0, 0);
    Point q2(0.2, 0.2);
    Point q3(0.4, 0.4);
    Point q4(0.6, 0.6);
    Point q5(0.8, 0.8);
    Point q6(1, 1);

    Curve c1("", {
        p1, p2, p3, p4, p5, p6
    });

    Curve c2("", {
        q1, q2, q3, q4, q5, q6
    });

    using _CDTW = CDTW<2, Norm::L1, Norm::L1>;
    auto cdtw = _CDTW(c1, c2);

    std::cout << "test case 4 cost: " << cdtw.cost() << std::endl;
}

void test_case_5() {
    Point p1(0, 0);
    Point p2(0, 0.1);
    Point p3(0, 0.2);
    Point p4(0, 0.3);
    Point p5(0, 0.4);
    Point p6(0, 0.5);
    Point p7(0, 0.6);
    Point p8(0, 0.7);
    Point p9(0, 0.8);
    Point p10(0, 0.9);
    Point p11(0, 1);

    Point q1(0, 0);
    Point q2(0.1, 0.1);
    Point q3(0.2, 0.2);
    Point q4(0.3, 0.3);
    Point q5(0.4, 0.4);
    Point q6(0.5, 0.5);
    Point q7(0.6, 0.6);
    Point q8(0.7, 0.7);
    Point q9(0.8, 0.8);
    Point q10(0.9, 0.9);
    Point q11(1, 1);

    Curve c1("", {
        p1, p2, p3, p4, p5, p6,
        p7, p8, p9, p10, p11
    });

    Curve c2("", {
        q1, q2, q3, q4, q5, q6,
        q7, q8, q9, q10, q11
    });

    using _CDTW = CDTW<2, Norm::L1, Norm::L1>;
    auto cdtw = _CDTW(c1, c2);

    std::cout << "test case 5 cost: " << cdtw.cost() << std::endl;
}

void test_case_6() {
    Point p1(0, 0);
    Point p2(0, 1);

    Point q1(0, 0);
    Point q2(1, 1);


    auto c1 = Curve("", {
        p1, p2
    });

    auto c2 = Curve("", {
        q1, q2
    });

    using _CDTW = CDTW<2, Norm::L1, Norm::L1>;
    auto cdtw = _CDTW(c1, c2);

    std::cout << "test case 6 cost: " << cdtw.cost() << std::endl;
}

void test_case_7() {
    Point p1(0, 0);
    Point p2(0, 1);

    Point q1(0, 1);
    Point q2(0, 0);

    Curve c1("", {p1, p2});
    Curve c2("", {q1, q2});

    using _CDTW = CDTW<2, Norm::L1, Norm::L1>;
    auto cdtw = _CDTW(c1, c2);

    std::cout << "test case 7 cost: " << cdtw.cost() << std::endl;
}

void test_case_8() {
    Point p1(0, 0);
    // Point p2(0, .2);
    // Point p3(0, .6);
    Point p4(0, .8);
    Point p5(0, .9);
    Point p6(0, 1);

    Point q1(0, 0);
    // Point q2(.2, .2);
    // Point q3(.6, .6);
    Point q4(.8, .8);
    Point q5(.9, .9);
    Point q6(1, 1);

    Curve c1("", {p1, /*p2, p3,*/ p4, p5, p6});
    Curve c2("", {q1, /*q2, q3,*/ q4, q5, q6});

    using _CDTW = CDTW<2, Norm::L1, Norm::L1>;
    auto cdtw = _CDTW(c1, c2);

    std::cout << "test case 8 cost: " << cdtw.cost() << std::endl;
}

void test_case_9() {
    Point p1(0, 0);
    Point p2(0, 0.5);
    Point p3(0, 1);

    Point q1(0, 0);
    Point q2(0.9, 0.9);
    Point q3(1, 1);

    Curve c1("", {p1, p2, p3});
    Curve c2("", {q1, q2, q3});

    using _CDTW = CDTW<2, Norm::L1, Norm::L1>;
    auto cdtw = _CDTW(c1, c2);

    std::cout << "test case 9 cost: " << cdtw.cost() << std::endl;
}

void test_case_10() {

    Point p1(0, 0);
    Point p2(0, 1);
    Point p3(0, 2);

    Point q1(0, 0);
    Point q2(0, 1);
    Point q3(0, 4);


    Curve c1("", {p1, p2, p3});
    Curve c2("", {q1, q2, q3});

    using _CDTW = CDTW<2, Norm::L1, Norm::L1>;
    auto cdtw = _CDTW(c1, c2);

    std::cout << "test case 10 cost: " << cdtw.cost() << std::endl;

}


int main() {
    // TODO: Compare Heuristic CDTW vs Exact CDTW (timing and result)
    // TODO: Upgrade to 2D + L2^2
    // TODO: How to Dijkstra in Exact CDTW?
    // auto curve1 = io::read_curve("data/characters/data/a0001.txt");
    // auto curve2 = io::read_curve("data/characters/data/z2796.txt");

    // auto curves = io::read_curves("data/characters/data");

    // std::cout << "number of curves: " << curves.size() << std::endl;

    Point p1(0, 0);
    Point p2(0, 1);
    Point p3(0, 0);
    Point p4(1, 1);


    // Curve c1("", {p1, p2});
    // Curve c2("", {p3, p4});

    // auto curve1_ = curve1.slice(1, 30);
    // auto curve2_ = curve2.slice(1, 30);
    
    // using _CDTW = CDTW<2, Norm::L1, Norm::L1>;
    // auto cdtw = _CDTW(c1, c2);

    // for (int i = 0; i <= 10; ++i) {
        // for (int j = 500; j <= 510; ++j) {
            // auto c_1 = curves[i].slice(1, 30);
            // auto c_2 = curves[j].slice(1, 30);
            // auto _cdtw = _CDTW(c_1, c_2);
            // std::cout << "cost: " << _cdtw.cost() << std::endl;
        // }
    // }

    // IntegralFrechet heuristic_alg(curve1_, curve2_, ParamMetric::L1, .1);
    //     const auto heuristic_res = heuristic_alg.compute_matching();

    // std::cout << heuristic_res.cost << std::endl;
    // auto cdtw = _CDTW(c1, c2);
    // std::cout << cdtw.cost() << std::endl;

    // test_case_1();
    // test_case_2();
    // test_case_3();
    // test_case_4();
    // test_case_5();
    // test_case_6();
    // test_case_7();
    // test_case_8();
    test_case_9();
    // test_case_10();

    Polynomial<1> left = Polynomial<1>(
            {{0, 0}}
    );
    Polynomial right = Polynomial<1>(
        {{.1, 0}}
    );
    Polynomial right_y = Polynomial<1>(
        {{0, 1}}
    );

    ConstrainedBivariatePolynomial<2> au = ConstrainedBivariatePolynomial<2>();
    au.f = BivariatePolynomial<2>({{
	{0.14142135623730948, 1.0999999999999999, -0.5},
	{-1.2414213562373093, 0, 0},
	{1.9142135623730951, 0, 0}
    }});
    
    au.path_type = AU;
    au.y_interval.max = .1;
    au.y_interval.min = 0;
    au.left_constraints.push_back(left);
    au.right_constraints.push_back(right);
    au.right_constraints.push_back(right_y);

    ConstrainedBivariatePolynomial<2> uau = ConstrainedBivariatePolynomial<2>();
    uau.f = BivariatePolynomial<2>({{
	    {{0.14142135623730948, 0.75857864376269046, 1.2071067811865475}},
	    {{-.9, 0, 0}},
	    {{0.20710678118654768, 0, 0}}
    }});

    uau.path_type = UAU;
    uau.y_interval.min = 0;
    uau.y_interval.max = .1;
    uau.left_constraints.push_back(left);
    uau.right_constraints.push_back(right);
    uau.right_constraints.push_back(right_y);

    for (int i = 0; i < 100; ++i) {
        for (int j = i+1; j <= 100; ++j) {
            double x = 0.001*i;
            double y = 0.001*j;

            double diff = au.slice_at_y(y).polynomial(x) - uau.slice_at_y(y).polynomial(x);
            if (diff < 0)
                std::cout << "that's a problem...\n";

        }
    }

    // std::cout << "au(0, .05) = " << au.slice_at_y(.05).polynomial(0) << std::endl;
    // std::cout << "uau(0,.05) = " << uau.slice_at_y(.05).polynomial(0) << std::endl;


    
    std::vector<ConstrainedBivariatePolynomial<2>> costs = {au, uau};

    auto in_p = PolynomialPiece<2>({0, .1}, Polynomial<2>({{0, 0, 0}}));

    PiecewisePolynomial<2> in = PiecewisePolynomial<2>();

    in.pieces.push_back(in_p);

    // auto polynomials = propagate(costs, in.pieces.begin(), in.pieces.end());

    // experiment_compare_heuristic_vs_extact_cdtw();

    return 0;
}
