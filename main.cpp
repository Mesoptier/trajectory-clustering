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

    void simplification_comparison_and_plot() {
        std::cout << std::endl;
        std::cout << std::string(80, '-') << "\nSIMPLIFICATION EXPERIMENT\n"
            << std::string(80, '-') << "\n";
        std::cout << "Load dataset, sample curves, subsample to complexity 50.\n"
            << "Run different simplification approaches.\n"
            << "For each curve, evaluate with CDTW and compute statistics over "
            << "the dataset.\n"
            << "Time is total time in ms. The rest are CDTW values.\n"
            << std::endl
            << "Save the simplifications and subsampled curves in ./simpl/.\n"
            << "Report the costs for the specific curves.\n"
            << "The curves chosen are the worst ones for specific approaches.\n"
            << std::endl;
        std::cout << "Character dataset\n";
        auto const curves_char = io::read_curves("data/characters/data");
        auto const samples_char = experiments::sample(curves_char, 10, 50);
        experiments::evaluate(samples_char, 12, "char");

        std::cout << "Pigeon dataset\n";
        auto const curves_pigeon = io::read_curves("data/Data_for_Mann_et_al_RSBL"
            "/Bladon & Church route recapping/bladon heath/utm", 1);
        auto const samples_pigeon = experiments::sample(curves_pigeon, 1, 40);
        experiments::evaluate(samples_pigeon, 6, "pigeon", 2);
    }

    void main_clustering_experiments() {
        std::cout << std::endl;
        std::cout << std::string(80, '-') << "\nINITIAL CLUSTERING EXPERIMENT\n"
            << std::string(80, '-') << "\n";
        std::cout << "Load dataset, sample curves, subsample to complexity 200.\n"
            << "Run different clustering approaches.\n"
            << "For each cluster, evaluate with k-medians + CDTW.\n"
            << "Save results in results/initial_clustering_experiment/"
            << "results.txt\n" << std::endl;
        experiments::initial_clustering_experiment();

        std::cout << std::endl;
        std::cout << std::string(80, '-') << "\nCENTER UPDATE EXPERIMENTS\n"
            << std::string(80, '-') << "\n";
        std::cout << "For each of the datasets, subsample curves to complexity 200.\n"
            << "After clustering with PAM, run different center update approaches.\n"
            << "For each cluster, evaluate with k-medians + CDTW.\n"
            << "Save results in results/ subdirectories.\n" << std::endl;
        std::cout << "Character dataset\n";
        // experiments::center_update_experiment_characters("characters", 50, 2, 8);
        std::cout << "Pigeon dataset\n";
        experiments::center_update_experiment_pigeons("pigeons", 3, 10);
        std::cout << "Movebank dataset\n";
        experiments::center_update_experiment_movebank("movebank", 3, 14);
    }

    void k_clustering_experiments() {
        std::cout << std::endl;
        std::cout << std::string(80, '-')
            << "\nCENTER UPDATE EXPERIMENTS FOR RANGE OF K\n"
            << std::string(80, '-') << "\n";
        std::cout << "For each dataset, do the center update experiments for "
            << "different values of k.\n"
            << "Save results in results/ subdirectories.\n" << std::endl;
        std::cout << "Character dataset\n";
        experiments::curve_complexity_experiment_characters();
        std::cout << "Pigeon dataset\n";
        experiments::curve_complexity_experiment_pigeons();
    }

    void extra_experiments() {
        std::cout << std::endl;
        std::cout << std::string(80, '-') << "\nSYNTHETIC CURVES EXPERIMENT\n"
            << std::string(80, '-') << "\n";
        std::cout << "Load the first handwritten a, subsample to complexity 200.\n"
            << "Perturb the points a bit to create many different curves.\n"
            << "Attempt to recreate the true center with a clustering approach.\n"
            << "Save results in plot.txt\n" << std::endl;
        experiments::synthetic_curve_experiment();

        std::cout << std::endl;
        std::cout << std::string(80, '-') << "\nCLASSIFICATION EXPERIMENT\n"
            << std::string(80, '-') << "\n";
        std::cout << "Load a sample of the character dataset.\n"
            << "Cluster the curves to check if we correctly identify the letters.\n"
            << "Report results on stdout.\n" << std::endl;
        classification::characterClassification();

        std::cout << std::endl;
        std::cout << std::string(80, '-') << "\nWEDGE PARAMETER EXPERIMENT\n"
            << std::string(80, '-') << "\n";
        std::cout << "Load the trajectories of a specific pigeon.\n"
            << "Identify the best parameters (eps, r) for the grid search.\n"
            << "Report results on stdout.\n" << std::endl;
        experiments::find_wedge_params_pigeons();
    }
}

int main() {
    thesis_example();
    // warping_path_example();
    // correctness_experiments();
    // negative_valley_test();
    // generate_figure();
    // running_time::experiment();
    // test_for_profiler();
    exit(0);
    static_assert(std::numeric_limits<distance_t>::is_iec559,
        "IEEE 754 arithmetic is assumed.");
    std::cout << "This is the codebase illustrating the experiments in the "
        << "paper\n(k, l)-Medians Clustering of Trajectories Using Continuous "
        << "Dynamic Time Warping\nby Milutin Brankovic, Kevin Buchin, Koen "
        << "Klaren, Andre Nusser, Aleksandr Popov,\nand Sampson Wong.\nPlease "
        << "consult https://github.com/Mesoptier/trajectory-clustering for more"
        << " info.\n\nWe first run the simplification-related experiments (see "
        << "Sections 5.2, 6.3);\nthen we evaluate the initial clustering "
        << "approaches on the pigeon dataset (see\nSections 5.1, 6.2); then we "
        << "compare the center update methods on all three\ndatasets, taking "
        << "the best simplification and initial clustering approaches (see\n"
        << "Sections 5.3, 6.4). Finally, we run some extra experiments that "
        << "may provide some\nextra information, but are less interesting and "
        << "did not make the cut.\n" << std::endl;
    // simplification_comparison_and_plot();
    main_clustering_experiments();
    char answer = 'z';
    do {
        std::cout << "Would you like to perform the center update experiments "
            << "on all datasets for\ndifferent values of k? This might take a "
            << "while. [y/n]\n";
        std::cin >> answer;
    } while (!std::cin.fail() && answer != 'y' && answer != 'n');
    if (answer == 'y')
        k_clustering_experiments();

    answer = 'z';
    do {
        std::cout << "Would you like to perform extra experiments (not described"
            << " in the paper)?\nThey are probably less interesting / incomplete."
            << " [y/n]\n";
        std::cin >> answer;
    } while (!std::cin.fail() && answer != 'y' && answer != 'n');
    if (answer == 'y')
        extra_experiments();
    std::cout << "Bye!" << std::endl;
    return 0;
}
