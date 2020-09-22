#include <iostream>
#include <string>

#include "src/Curve.h"
#include "src/classification_experiment.h"
#include "src/experiments.h"
#include "src/simplification_experiment.h"
#include "src/utils/io.h"

namespace {
    void simplification_comparison_and_plot() {
        std::cout << "\n" << std::endl;
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
        auto const curves_pigeon = io::read_curves("data/Data_for_Mann_et_al_RSBL 2"
            "/Bladon & Church route recapping/bladon heath/utm", 1);
        auto const samples_pigeon = experiments::sample(curves_pigeon, 1, 40);
        experiments::evaluate(samples_pigeon, 6, "pigeon", 2);
    }

    void main_clustering_experiments() {
        std::cout << "\n" << std::endl;
        std::cout << std::string(80, '-') << "\nINITIAL CLUSTERING EXPERIMENT\n"
            << std::string(80, '-') << "\n";
        std::cout << "Load dataset, sample curves, subsample to complexity 200.\n"
            << "Run different clustering approaches.\n"
            << "For each cluster, evaluate with k-medians + CDTW.\n"
            << "Save results in results/initial_clustering_experiment/"
            << "results.txt\n" << std::endl;
        experiments::initial_clustering_experiment();

        std::cout << "\n" << std::endl;
        std::cout << std::string(80, '-') << "\nCENTER UPDATE EXPERIMENTS\n"
            << std::string(80, '-') << "\n";
        std::cout << "For each of the datasets, subsample curves to complexity 200.\n"
            << "After clustering with PAM, run different center update approaches.\n"
            << "For each cluster, evaluate with k-medians + CDTW.\n"
            << "Save results in results/ subdirectories.\n" << std::endl;
        std::cout << "Character dataset\n";
        experiments::center_update_experiment_characters("characters", 50, 2, 8);
        std::cout << "Pigeon dataset\n";
        experiments::center_update_experiment_pigeons("pigeons", 3, 10);
        std::cout << "Movebank dataset\n";
        experiments::center_update_experiment_movebank("movebank", 3, 14);
    }

    void k_clustering_experiments() {
        std::cout << "\n" << std::endl;
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
        std::cout << "\n" << std::endl;
        std::cout << std::string(80, '-') << "\nSYNTHETIC CURVES EXPERIMENT\n"
            << std::string(80, '-') << "\n";
        std::cout << "Load the first handwritten a, subsample to complexity 200.\n"
            << "Perturb the points a bit to create many different curves.\n"
            << "Attempt to recreate the true center with a clustering approach.\n"
            << "Save results in plot.txt\n" << std::endl;
        experiments::synthetic_curve_experiment();

        std::cout << "\n" << std::endl;
        std::cout << std::string(80, '-') << "\nCLASSIFICATION EXPERIMENT\n"
            << std::string(80, '-') << "\n";
        std::cout << "Load a sample of the character dataset.\n"
            << "Cluster the curves to check if we correctly identify the letters.\n"
            << "Report results on stdout.\n" << std::endl;
        classification::characterClassification();

        std::cout << "\n" << std::endl;
        std::cout << std::string(80, '-') << "\nWEDGE PARAMETER EXPERIMENT\n"
            << std::string(80, '-') << "\n";
        std::cout << "Load the trajectories of a specific pigeon.\n"
            << "Identify the best parameters (eps, r) for the grid search.\n"
            << "Report results on stdout.\n" << std::endl;
        experiments::find_wedge_params_pigeons();
    }
}

int main() {
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
    simplification_comparison_and_plot();
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
