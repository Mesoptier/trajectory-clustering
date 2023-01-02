#include <iostream>
#include <string>
#include <stdlib.h>

#include "src/Curve.h"
#include "src/classification_experiment.h"
#include "src/experiments.h"
#include "src/simplification_experiment.h"
#include "src/utils/io.h"
#include "src/cdtw/2d-l1-l1.h"
#include "src/cdtw/cdtw.h"
#include "src/test.h"
#include "src/running_time_exp.h"
#include <random>


namespace {
    void value_test() {
        auto const ch_curves = io::read_curves("data/characters/data");
        auto const curves = experiments::sample(ch_curves, 10, 50);

        auto c1 = curves[1].naive_l_simplification(100);
        auto c2 = curves[200].naive_l_simplification(100);

        auto start_time = std::chrono::high_resolution_clock::now();
        auto val = df::cdtw_2d_l1_l1(c1, c2);
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
        std::cout << "time: " << duration << std::endl;
        std::cout << val << std::endl; // should be 14.0511
    }

    void value_test_2() {
        Points p1 {
            {0, 0}, {1, 1}, {2, 0}
        };

        Points p2 {
            {0, 1}, {1, 2}, {2, 1}
        };

        Curve c1{"", p1};
        Curve c2{"", p2};

        auto cdtw = CDTW<2, Norm::L1, Norm::L1>(c1, c2);
        std::cout << cdtw.cost() << "\n";
    }

    void value_test_3() {
        Points p1 {
            {0, 0}, {1, 1}, {2, 0}
        };

        Points p2 {
            {0, 1}, {0, 2}, {0, 3}
        };

        Curve c1{"", p1};
        Curve c2{"", p2};

        auto cdtw = CDTW<2, Norm::L1, Norm::L1>(c1, c2);
        std::cout << cdtw.cost() << "\n";
    }

    void value_test_4() {
        Points p1 {
            {0, 0}, {2, 0}, {3, 0}
        };

        Points p2 {
            {0, 1}, {0, 2}, {0, 3}
        };

        Curve c1{"", p1};
        Curve c2{"", p2};

        auto cdtw = CDTW<2, Norm::L1, Norm::L1>(c1, c2);
        std::cout << cdtw.cost() << "\n";
    }

    void value_test_5() {
        Points p1 {
            {0, 0}, {1, 1}
        };

        Points p2 {
            {0, 1}, {1, 0}
        };

        Curve c1{"", p1};
        Curve c2{"", p2};

        auto cdtw = CDTW<2, Norm::L1, Norm::L1>(c1, c2);
        std::cout << cdtw.cost() << "\n";
    }


    void thesis_example() {
        Points p1 {
            {0, 0}, {0, 2}
        };

        Points p2 {
            {1,1}, {2,1}
        };
        
        Curve c1{"", p1};
        Curve c2{"", p2};

        auto cdtw = CDTW<2, Norm::L1, Norm::L1>(c1, c2);

        std::cout << cdtw.cost() << "\n";
    }


    Curves get_rand_order_simplified_curves(int res, std::string dataset) {
        auto const pi_curves = dataset == "pigeons" ? io::read_curves("data/Data_for_Mann_et_al_RSBL"
            "/Bladon & Church route recapping/bladon heath/utm", 1) : 
            io::read_curves("data/characters/data");

        std::vector<int> pi_indices = std::vector<int>();

        for (int i = 0; i < pi_curves.size(); ++i) 
            pi_indices.push_back(i);

        int seed = 5;
        std::mt19937_64 g(seed);

        std::shuffle(pi_indices.begin(), pi_indices.end(), g);

        Curves curves = Curves();

        if (res > -1)
            for (auto i: pi_indices)
                curves.push_back(pi_curves[i].naive_l_simplification(res));
        else
            for (auto i: pi_indices)
                curves.push_back(pi_curves[i]);

        return curves;
    }

    void warping_path_example () {
        Curves curves = get_rand_order_simplified_curves(-1, "character");

        auto heur_cost = df::heur_cdtw_2d_l1_l1_path(curves[0], curves[1]);
        auto cdtw = CDTW<2, Norm::L1, Norm::L1>(curves[0], curves[1]);
        auto cost = cdtw.cost();
        cdtw.write_warping_path();

        std::cout << "heur_cost: " << heur_cost << std::endl;
        std::cout << "exact cost: " << cost << std::endl;
    }

    void correctness_experiments() {

        auto const ch_curves = io::read_curves("data/characters/data");

        auto const pi_curves = io::read_curves("data/Data_for_Mann_et_al_RSBL"
            "/Bladon & Church route recapping/bladon heath/utm", 1);

        std::vector<int> ch_indices = std::vector<int>();

        for (int i = 0; i <= ch_curves.size(); ++i)
            ch_indices.push_back(i);

        std::vector<int> pi_indices = std::vector<int>();

        for (int i = 0; i <= pi_curves.size(); ++i) 
            pi_indices.push_back(i);

        int seed = 10;
        std::mt19937_64 g(seed);
        std::shuffle(ch_indices.begin(), ch_indices.end(), g);
        std::shuffle(pi_indices.begin(), pi_indices.end(), g);

        std::ofstream ch_results("2d_paper_results/correctness/character_results.csv");
        std::ofstream pi_results("2d_paper_results/correctness/pigeon_results.csv");

        ch_results << "exact,heuristic\n";
        pi_results << "exact,heuristic\n";


        for (int i = 1; i <= 39; i+=2) {
            std::cout << i << std::endl;

            Curve c1 = ch_curves[ch_indices[i]].naive_l_simplification(100);
            Curve c2 = ch_curves[ch_indices[i+1]].naive_l_simplification(100);

            auto heur_cost = df::heur_cdtw_2d_l1_l1(c1, c2);
            auto exact_cost = CDTW<2, Norm::L1, Norm::L1>(c1, c2).cost();
            ch_results << exact_cost << "," << heur_cost << std::endl;
        }

        for (int i = 1; i <= 39; i+=2) {
            std::cout << i << std::endl;

            Curve c1 = pi_curves[pi_indices[i]].naive_l_simplification(100);
            Curve c2 = pi_curves[pi_indices[i+1]].naive_l_simplification(100);

            auto heur_cost = df::heur_cdtw_2d_l1_l1(c1, c2);
            auto exact_cost = CDTW<2, Norm::L1, Norm::L1>(c1, c2).cost();
            pi_results << exact_cost << "," << heur_cost << std::endl;
        }

        ch_results.close();
        pi_results.close();
    }

    void negative_valley_test() {

        Points p1 = {
            {.55, 0},
            {1, 0}
        };

        Points p2 = {
            {1, .2},
            {0, -.2}
        };

        Curve c1("", p1);
        Curve c2("", p2);

        auto cdtw = CDTW<2, Norm::L1, Norm::L1>(c1, c2);
        std::cout << "cdtw cost:" << cdtw.cost() << std::endl;

        std::cout << "heursitic cost: " << df::heur_cdtw_2d_l1_l1(c1, c2) << std::endl;
    }

    void generate_figure() {

        Points p1 = {
            {1, 1},
            {2, 1.1},
            {3, 0.9},
        };

        Points p2 = {
            {1, 2},
            {2, 1.9},
            {2.5, 1.5},
            {3, 2},
            {3.5, 1.1}
        };

        Curve c1("", p1);
        Curve c2("", p2);

        auto cdtw = CDTW<2, Norm::L1, Norm::L1>(c1, c2);
        std::cout << cdtw.cost() << std::endl;
    }

    void test_for_profiler() {
        auto const ch_curves = io::read_curves("data/characters/data");
        auto const curves = experiments::sample(ch_curves, 10, 50);

        auto val = df::cdtw_2d_l1_l1(curves[1].naive_l_simplification(100), curves[200].naive_l_simplification(100));
        std::cout << val << std::endl;
    }

    void test_l1() {
        auto const curves = io::read_curves("data/characters/data");
        std::cout << "exact: " << df::cdtw_2d_l1_l1(curves[0].naive_l_simplification(10), curves[200].naive_l_simplification(10)) << std::endl;
        std::cout << "heur: " << df::heur_cdtw_2d_l1_l1(curves[0].naive_l_simplification(10), curves[200].naive_l_simplification(10)) << std::endl;
    }

    void lower_envelope_test() {
        std::vector<PolynomialPiece<2>> polynomials = std::vector<PolynomialPiece<2>>();
        int count = 1000;
        for (int i = 0; i < count; ++i) {
            double root1 = rand() % 10;
            // double root2 = rand() % 10;
            double root2 = root1 + .1;
            // double coeff1 = 1 + rand() % 9;
            // double coeff2 = 1 + rand() % 9;
            double coeff1 = 1;
            double coeff2 = 1;
            // double sign = rand() % 2 < .5 ? -1. : 1.;
            double sign = 1;
            double scale = 1 + rand() % 1000000;
            double a = sign*(root1 * root2 / (coeff1 * coeff2));
            double b = sign*(-root2 / coeff1 - root1 / coeff2);
            double c = sign * coeff1 * coeff2;

            auto polynomial = Polynomial<2>({{a, b, c}});
            // std::cout << polynomial << std::endl;
            double start = rand() % 10;
            double length = 1 + rand() % 2;
            Interval_c domain = {start, start + length};
            auto piece = PolynomialPiece<2>(domain, polynomial);
            polynomials.push_back(piece);
        }

        // auto envelope = naive_lower_envelope(polynomials);
        auto envelope = fast_lower_envelope_v2(polynomials);

        write_polynomial_set(polynomials, "test_data/polynomial_set.txt");
        write_polynomial_set(envelope.pieces, "test_data/lower_envelope.txt");
    }

}


void debug() {
    std::string filename = "test_data/polynomial_set.txt";
    auto pieces = read_pieces<2>(filename);
    auto envelope = fast_lower_envelope_v2(pieces);
    write_polynomial_set(envelope.pieces, "test_data/lower_envelope.txt");
}

int main() {
    value_test();
    // value_test_2();
    // value_test_3();
    // value_test_4();
    // value_test_5();
    // thesis_example();
    // warping_path_example();
    // correctness_experiments();
    // negative_valley_test();
    // generate_figure();
    // running_time::experiment();
    // test_for_profiler();
    // debug();
    // lower_envelope_test();
    exit(0);
}
