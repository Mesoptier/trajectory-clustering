#include <iostream>
#include <string>

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

        auto val = df::cdtw_2d_l1_l1(curves[1].naive_l_simplification(100), curves[200].naive_l_simplification(100));
        std::cout << val << std::endl;

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

        auto heur_cost = df::heur_cdtw_2d_l1_l1(curves[0], curves[1]);
        auto cost = CDTW<2, Norm::L1, Norm::L1>(curves[0], curves[1]).cost();

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


        for (int i = 1; i <= 41; i+=2) {
            std::cout << i << std::endl;

            Curve c1 = ch_curves[ch_indices[i]].naive_l_simplification(100);
            Curve c2 = ch_curves[ch_indices[i+1]].naive_l_simplification(100);

            auto heur_cost = df::heur_cdtw_2d_l1_l1(c1, c2);
            auto exact_cost = CDTW<2, Norm::L1, Norm::L1>(c1, c2).cost();
            ch_results << exact_cost << "," << heur_cost << std::endl;
        }

        for (int i = 1; i <= 41; i+=2) {
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

}

int main() {
    // value_test();
    // thesis_example();
    warping_path_example();
    // correctness_experiments();
    // negative_valley_test();
    // generate_figure();
    // running_time::experiment();
    // test_for_profiler();
    exit(0);
}
