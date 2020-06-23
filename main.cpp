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
#include "src/greedy_simplification.h"
#include "src/curve_simplification.h"
#include "src/CurveSimpMatrix.h"
#include "src/clustering/pam_with_simplifications.h"
#include "src/clustering/clustering_algs.h"
#include "src/clustering/center_algs.h"
#include "src/greedy_l_simplification.h"
#include "src/simplification/imaiiri.h"
#include "src/DTW/dtw.h"
#include "src/distance_functions.h"
#include "src/experiments.h"
#include "src/synthetic_curves.h"
#include "src/classification_experiment.h"
 


namespace {
//
// I/O helpers
//

std::vector<Curve> read_curves(const std::string& filename) {
    auto start_time = std::chrono::high_resolution_clock::now();

    const auto curves = io::read_curves(filename);

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    std::cout << "read_curves: " << duration << "ms\n";

    return curves;
}

SymmetricMatrix read_matrix(const std::string& filename) {
    auto start_time = std::chrono::high_resolution_clock::now();

    std::ifstream file(filename);
    const auto matrix = SymmetricMatrix::read(file);

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    std::cout << "read_matrix: " << duration << "ms\n";

    return matrix;
}

void export_matrix(const SymmetricMatrix& matrix, const std::string& filename) {
    std::ofstream file(filename);
    matrix.write(file);
    file.close();
}


//
// Algorithms
//

SymmetricMatrix compute_distance_matrix(const std::vector<Curve>& curves) {
    size_t n = curves.size();
    SymmetricMatrix distance_matrix(n);

    auto start_time = std::chrono::high_resolution_clock::now();

    for (size_t i = 0; i < n; ++i) {
        const auto curve1 = curves[i];

        // Don't compute the j == i case, since the cost is always 0 (comparing identical curves)
        distance_matrix.at(i, i) = 0;

        for (size_t j = i + 1; j < n; ++j) {
            const auto curve2 = curves[j];

            IntegralFrechet alg(curve1, curve2, ParamMetric::LInfinity_NoShortcuts, 100, nullptr);
            const auto result = alg.compute_matching();

            // TODO: Weigh by max or sum?
            // TODO: Move to IntegralFrechet class
            distance_t averaged_cost = result.cost / std::max(curve1.curve_length(), curve2.curve_length());
            distance_matrix.at(i, j) = averaged_cost;
        }

        std::cout << i << '\n';
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    std::cout << "compute_distance_matrix: " << duration << "ms\n";

    return distance_matrix;
}

void compute_clusters(const SymmetricMatrix& distance_matrix, size_t k) {
    auto start_time = std::chrono::high_resolution_clock::now();

    clustering::pam::compute(distance_matrix.n, k, distance_matrix);

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    std::cout << "compute_clusters: " << duration << "ms\n";
}

IntegralFrechet::MatchingResult compute_matching(const Curve& curve1, const Curve& curve2, distance_t resolution, const MatchingBand* const band) {
    auto start_time = std::chrono::high_resolution_clock::now();

    IntegralFrechet alg(curve1, curve2, ParamMetric::LInfinity_NoShortcuts, resolution, band);
    const auto result = alg.compute_matching();

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    std::cout << "compute_matching: " << duration << "ms\n";

    return result;
}

MatchingBand compute_band(const Curve& curve1, const Curve& curve2, const Points& matching, distance_t radius) {
    auto start_time = std::chrono::high_resolution_clock::now();

    const MatchingBand band(curve1, curve2, matching, radius);

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    std::cout << "compute_band: " << duration << "ms\n";

    return band;
}

//
// Experiments
//

void experiment_with_or_without_bands() {
    struct Stat {
        distance_t band_cost;
        long long int band_time;
        distance_t no_band_cost;
        long long int no_band_time;
    };
    std::vector<Stat> stats;

    bool maintain_lengths = true;
    distance_t simple_resolution = 10;
    distance_t resolution = 1;
    distance_t band_radius = 1;

    const auto curves = read_curves("data/characters/data");

    std::random_device rd;
    std::mt19937_64 generator(rd());
    std::uniform_int_distribution<size_t> distribution(0, curves.size() - 1);

    size_t k_max = 1000;
    for (size_t k = 0; k < k_max; ++k) {
        auto i = distribution(generator);
        auto j = distribution(generator);

        const auto curve1 = curves.at(i);
        const auto curve2 = curves.at(j);

        // With band
        auto start_band = std::chrono::high_resolution_clock::now();
        const auto result_alt = compute_matching(curve1.simplify(maintain_lengths), curve2.simplify(maintain_lengths), simple_resolution,nullptr);
        const auto band = compute_band(curve1, curve2, result_alt.matching, band_radius);
        const auto result_band = compute_matching(curve1, curve2, resolution, &band);
        auto end_band = std::chrono::high_resolution_clock::now();
        auto duration_band = std::chrono::duration_cast<std::chrono::milliseconds>(end_band - start_band).count();

        // Without band
        auto start_no_band = std::chrono::high_resolution_clock::now();
        const auto result_no_band = compute_matching(curve1, curve2, resolution, nullptr);
        auto end_no_band = std::chrono::high_resolution_clock::now();
        auto duration_no_band = std::chrono::duration_cast<std::chrono::milliseconds>(end_no_band - start_no_band).count();

        stats.push_back({
            result_band.cost,
            duration_band,
            result_no_band.cost,
            duration_no_band,
        });

        std::cout << (k + 1) << " / " << k_max << '\n';
    }

    std::ofstream file("data/out/band_stats.tsv");
    file << "band_cost\tband_time\tno_band_cost\tno_band_time\n";
    for (const auto& stat : stats) {
        file << stat.band_cost << '\t' << stat.band_time << '\t' << stat.no_band_cost << '\t' << stat.no_band_time << '\n';
    }
    file.close();
}

void experiment_visualize_band() {
    bool maintain_lengths = true;

    const auto curve1 = io::read_curve("data/characters/data/a0001.txt");
    const auto curve2 = io::read_curve("data/characters/data/a0002.txt");

    const auto result_no_band = compute_matching(curve1, curve2, 1, nullptr);
    io::export_points("data/out/debug_points_old.csv", result_no_band.search_stat.nodes_as_points);
    io::export_points("data/out/matching3.csv", result_no_band.matching);
    std::cout << "matching (no band) cost: " << result_no_band.cost << '\n';

    auto start_time = std::chrono::high_resolution_clock::now();

    const auto result_alt = compute_matching(curve1.simplify(maintain_lengths), curve2.simplify(maintain_lengths), 10,nullptr);
    const auto band = compute_band(curve1, curve2, result_alt.matching,     2);
    const auto result = compute_matching(curve1, curve2, 1, &band);

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    std::cout << "main: " << duration << "ms\n";

    io::export_points("data/out/curve1.csv", curve1.get_points());
    io::export_points("data/out/curve2.csv", curve2.get_points());
    io::export_points("data/out/matching1.csv", result.matching);
    io::export_points("data/out/matching2.csv", result_alt.matching);
    io::export_points("data/out/debug_points.csv", result.search_stat.nodes_as_points);

    std::cout << "matching cost: " << result.cost << '\n';
    std::cout << "matching (alt) cost: " << result_alt.cost << '\n';
}

std::vector<Curve> sample_curves(std::vector<Curve> curves, unsigned period) {

    std::vector<Curve> samples = std::vector<Curve>();

    for (std::size_t i = 0; i < curves.size(); ++i) {
        if (i % period == 0) {
            samples.push_back(curves[i]);
        }
    }

    return samples;
}

void evaluate_greedy_simplification() {
    std::vector<Curve> curves = sample_curves(read_curves("data/characters/data"), 20);

    std::vector<std::size_t> greedy_node_count;
    std::vector<distance_t> greedy_integral_distance;
    std::vector<std::size_t> benchmark_node_count;
    std::vector<distance_t> benchmark_integral_distance;

    std::fstream greedy;
    greedy.open("greedy_stats.txt", std::fstream::out | std::fstream::trunc);
    std::fstream benchmark;
    benchmark.open("benchmark_stats.txt", std::fstream::out | std::fstream::trunc);

    for (Curve curve: curves) {

        std::cout << curve.name() << "\n";

        Curve simplification = curve.simplify(false);
        Curve greedy_simp = greedy_simplification(curve, 0.25);

        // io::export_points(curve.name() + "-simplification.txt", greedy_simp.get_points());
        // io::export_points(curve.name() + "-regular-simplification.txt", simplification.get_points());

        greedy_node_count.push_back(greedy_simp.get_points().size());
        //std::cout << "added greedy count \n";
        greedy_integral_distance.push_back(
            IntegralFrechet(greedy_simp, curve, ParamMetric::LInfinity_NoShortcuts, 1, nullptr).compute_matching().cost
        );
        //std::cout << "added greedy distance \n";
        benchmark_node_count.push_back(simplification.get_points().size());
        //std::cout << "added benchmark count \n";
        benchmark_integral_distance.push_back(
            IntegralFrechet(simplification, curve, ParamMetric::LInfinity_NoShortcuts, 1, nullptr).compute_matching().cost
        );

        greedy << greedy_node_count.back() << "," << greedy_integral_distance.back() << "\n";
        benchmark << benchmark_node_count.back() << "," << benchmark_integral_distance.back() << "\n";
    }

    greedy.close();
    benchmark.close();
}



void evaluate_frechet_simplifications() {
    std::vector<Curve> curves = sample_curves(read_curves("data/characters/data"), 20);
    const auto simplifications = read_curves("data/simplifications");

    std::fstream frechet_stats;
    frechet_stats.open("frechet_stats.csv", std::fstream::out | std::fstream::trunc);

    for (std::size_t i = 0; i < curves.size(); ++i) {

        Curve original_curve = curves[i];
        Curve simplification = simplifications[i];

        const auto distance =
        IntegralFrechet(original_curve, simplification, ParamMetric::LInfinity_NoShortcuts, 1, nullptr)
        .compute_matching()
        .cost;

        frechet_stats << simplification.get_points().size() << "," << distance << "\n";
    }

    frechet_stats.close();
}

void write_simplifications() {

    Curves curves = read_curves("data/characters/data");

    Curves greedy_simplifications = Curves();
    Curves regular_simplifications = Curves();

    std::cout << "computing greedy simplifications... \n";
    for (auto curve: curves) {
        greedy_simplifications.push_back(
            greedy_simplification(curve, 0.25)
        );
    }

    std::cout << "computing regular simplifications... \n";
    for (auto curve: curves) {
        regular_simplifications.push_back(
            curve.simplify(true)
        );
    }

    std::fstream greedy_dataset_file;
    greedy_dataset_file.open("greedy_simplifications/dataset.txt", std::fstream::out | std::fstream::trunc);

    std::fstream regular_dataset_file;
    regular_dataset_file.open("regular_simplifications/dataset.txt", std::fstream::out | std::fstream::trunc);

    for (std::size_t i = 0; i < curves.size(); ++i) {

        greedy_dataset_file << "greedy_simplifications/gs-" + std::to_string(i) + ".txt" << "\n";        
        std::fstream greedy_simp;
        greedy_simp.open("greedy_simplifications/gs-" + std::to_string(i) + ".txt", std::fstream::out | std::fstream::trunc);
        for (auto point: greedy_simplifications[i].get_points()) {
            greedy_simp << point.x << " " << point.y << "\n";
        }
        greedy_simp.close();

        regular_dataset_file << "regular_simplifications/rs-" + std::to_string(i) + ".txt" << "\n";  
        std::fstream regular_simp;
        regular_simp.open("regular_simplifications/rs-" + std::to_string(i) + ".txt", std::fstream::out | std::fstream::trunc);
        for (auto point: regular_simplifications[i].get_points()) {
            regular_simp << point.x << " " << point.y << "\n";
        }
        regular_simp.close();
    }

    regular_dataset_file.close();
    greedy_dataset_file.close();

}


void test_clustering_algs() {
    // Curves curves = sample_curves(io::read_pigeon_curves("data/Data_for_Mann_et_al_RSBL 2/Bladon & Church route recapping/bladon heath"), 30);
    Curves curves = sample_curves(read_curves("data/characters/data"), 25);
    std::cout << "read curves...\n";
    std::cout << curves.size() << "\n";

    Clustering gonzalez_clustering = runGonzalez(curves, 26, 5, integral_frechet, true, "");
    std::cout << "finshed gonzalez...\n";

    Clustering single_linkage_clustering = singleLinkage(curves, 26, 5, integral_frechet, true);
    std::cout << "finished single linkage...\n";

    Clustering complete_linkage_clustering = completeLinkage(curves, 26, 5, integral_frechet, true);
    std::cout << "finished complete linkage...\n";

    distance_t gonzalez_sum = 0;
    for (auto cluster: gonzalez_clustering) {
        CurveIDs ids = cluster.curve_ids;
        for (auto id: ids) {
            gonzalez_sum += IntegralFrechet(curves[id], cluster.center_curve, ParamMetric::LInfinity_NoShortcuts, 1, nullptr)
            .compute_matching()
            .cost;
        }
    }

    std::cout << "gonzalez: " << gonzalez_sum << "\n";
    

    distance_t single_linkage_sum = 0;
    for (auto cluster: single_linkage_clustering) {
        for (auto id: cluster.curve_ids)
            single_linkage_sum += IntegralFrechet(curves[id], cluster.center_curve, ParamMetric::LInfinity_NoShortcuts, 1, nullptr)
            .compute_matching()
            .cost;
    }

    std::cout << "single linkage: " << single_linkage_sum << "\n";

    distance_t complete_linkage_sum = 0;
    for (auto cluster: complete_linkage_clustering) {
        for (auto id: cluster.curve_ids)
            complete_linkage_sum += IntegralFrechet(curves[id], cluster.center_curve, ParamMetric::LInfinity_NoShortcuts, 1, nullptr)
            .compute_matching()
            .cost;
    }

    std::cout << "complete linkage: " << complete_linkage_sum << "\n";
}

distance_t evaluate_clustering(Clustering clustering, Curves curves, distance_t(*dist_func)(Curve, Curve)) {

    distance_t cost = 0;

    for (auto cluster: clustering) {
        for (auto curve_id: cluster.curve_ids) {
            cost += dist_func(curves[curve_id], cluster.center_curve);
        }
    }

    return cost;
}

// void test_center_algs() {
//     Curves curves = read_curves("data/characters/data");
//     curves = Curves(curves.begin(), curves.begin() + 26);
//     /*
//     Curve curve = curves[0];
//     Curve segment = Curve("", {curve.get_points()[6], curve.get_points()[8]});
//     std::cout << curve.name() << "\n";
//     std::cout << IntegralFrechet(Curve("", {curve.get_points()[6], curve.get_points()[7], curve.get_points()[8]}), segment, ParamMetric::LInfinity_NoShortcuts, 100)
//     .compute_matching().cost << "\n";*/

//     Clustering clustering = runGonzalez(curves, 1, 10, integral_frechet, true);
//     std::cout << "computed initial cluster center...\n";
//     // clustering[0].center_curve = curves[0].simplify(true);

//     std::fstream script;
//     script.open("gnuplot_script.txt", std::fstream::out | std::fstream::trunc);

//     script << "plot ";
//     for (std::size_t i = 0; i < curves.size(); ++i) {
//         Curve curve = curves[i];
//         script << "\"" << curve.name() + "\" with linespoints ls 0.7 lt rgb \"black\" ps 0.01, ";
//     }

//     std::fstream initial_center;
//     initial_center.open("cluster/initial_center.txt", std::fstream::out | std::fstream::trunc);
//     for (auto p: clustering[0].center_curve.get_points()) {
//         initial_center << p.x << " " << p.y << "\n";
//     }
//     initial_center.close();

//     script << "\"cluster/initial_center.txt\" with linespoints ls 2 lw 3 lt rgb \"green\", ";

//     for (std::size_t i = 0; i < curves.size(); ++i) {
//         Curve curve = curves[i];
//         Points matching = matching_of_vertices(clustering[0].center_curve, curve);
//         std::fstream vertices;
//         vertices.open("cluster/matching" + std::to_string(i) + ".txt", std::fstream::out | std::fstream::trunc);
//         for (std::size_t j = 0; j < matching.size(); ++j) {
//             vertices << matching[j].x << " " << matching[j].y << "\n";
//             vertices << clustering[0].center_curve.get_points()[j].x << " " << clustering[0].center_curve.get_points()[j].y << "\n\n";
//         }
//         vertices.close();
//         script << "\"cluster/matching" + std::to_string(i) + ".txt\" with linespoints ls 1 lt rgb \"blue\", ";
//     }

//     int improvement_count = 0;
//     while (calcFSACenters(curves, clustering, 10, average_frechet, C2CDist::Median, CenterCurveUpdateMethod::frechetMean)) {
//         std::cout << "found new center!!\n";
//         improvement_count++;
//     }

//     std::cout << "finished\n";
//     std::cout << improvement_count << "\n";

//     std::fstream improved_cluster;
//     improved_cluster.open("cluster/improved_cluster.txt", std::fstream::out | std::fstream::trunc);
//     for (auto p: clustering[0].center_curve.get_points()) {
//         improved_cluster << p.x << " " << p.y << "\n";
//     }
//     improved_cluster.close();
//     script << "\"cluster/improved_cluster.txt\" with linespoints ls 3 lw 3 lt rgb \"red\"";
//     script.close();
// }

void test_pam_with_centering() {
    Curves curves = sample_curves(read_curves("data/characters/data"), 30);
    Clustering clustering = pam_with_centering(curves, 10, 10, integral_frechet, "");
}

void test_frechet() {
    Curves curves = read_curves("data/characters/data");
    curves = sample_curves(curves, 50);


    Curves simplifications = Curves();
    for (auto& curve: curves) {
        std::cout << curve.name() << "\n";
        simplifications.push_back(
            simplify(curve, 10, frechet)
        );
    }
    // Clustering pam_clustering = pam_with_centering(curves, 10, 10, average_frechet, "characters_matrix.txt");
    // Clustering gonzalez_clustering = runGonzalez(curves, 10, 10, frechet, true);

    // std::cout << "cost of pam clustering...\n";
    // std::cout << evaluate_clustering(pam_clustering, curves, average_frechet) << "\n";

    // std::cout << "cost of gonzalez clustering...\n";
    // std::cout << evaluate_clustering(gonzalez_clustering, curves, average_frechet) << "\n";
}

void compute_curve_simp_matrix() {
    Curves curves = read_data();
    Curves simplifications = Curves();

    for (auto curve: curves) {
        simplifications.push_back(
            curve.naive_l_simplification(10)
        );

        std::cout << simplifications.back().get_points().size() << "\n";
    }

    std::cout << "computed simplifications...\n";

    CurveSimpMatrix matrix = CurveSimpMatrix(curves, simplifications, average_frechet);
    matrix.write("pigeon_matrix.txt");
}

}

int main() {



//    const auto distance_matrix = read_matrix("data/out/distance_matrix.mtx");
//    compute_clusters(distance_matrix, 19); // "characters" has 19 classes

    //const auto curves = read_curves("data/characters/data");
    //const auto distance_matrix = compute_distance_matrix(curves);
    //export_matrix(distance_matrix, "data/out/distance_matrix.mtx");

    //experiment_visualize_band();
    //evaluate_greedy_simplification();
    // evaluate_frechet_simplifications();
    // test_clustering();
    // write_simplifications();
    // test_clustering_algs();
    // test_center_algs();
    // test_pam_with_centering();

    // const auto curves = read_curves("data/characters/data");
    // const auto dm = compute_distance_matrix(curves);
    // export_matrix(dm, "data/out/distance_matrix.mtx");
    // const auto distance_matrix = read_matrix("data/out/distance_matrix.mtx");
    // compute_clusters(distance_matrix, 19); // "characters" has 19 classes

    // std::vector<double> vec = {.001, .021, 1.12, 1.121, 1.23, 1.23456};
    // auto it = std::lower_bound(vec.begin(), vec.end(), 1.234);
    // auto index = static_cast<std::size_t>(std::distance(vec.begin(), it));
    // std::cout << index << "\n";



    curve_complexity_experiment_pigeons();
    // running_time_experiment();
    // compute_matrices();
    // full_clustering_experiment();
    // find_params();
    // pigeon_experiment();
    // wedge_method_experiment();
    // test_frechet();
    // compute_curve_simp_matrix();
    // run_experiments();
    // preliminary_experiments();
    // center_update_experiment_pigeons();
    // center_update_experiments();
    // ensemble_experiment();
    // characterClassification();
    // generate_curves(Curve({{0, 0}, {10, 10}, {20, 20}, {30, 30}, {40, 40}, {50, 50}, {60, 60}, {70, 70}, {80, 80}, {90, 90}, {100, 100}, {110, 110}, {120, 120}, {130, 130}}), 5);
    // write_curves();
    // synthetic_curve_experiment();
    // Curve curve1 = Curve({{-13.519655, 518.2176}, {-13.4540596224652, 518.171459590597}, {-13.4312466704293, 518.128483280568}, {-13.3730548740233, 518.11155638744}, {-13.3532487713532, 518.103879442256}, {-13.3012273722711, 518.081146702316}, {-13.2822606164468, 518.062349118219}, {-13.2407413511353, 518.010277368733}, {-13.2249619456808, 517.946611346145}, {-13.2343869456808, 517.908641346145}, {-13.2343869456808, 517.908641346145}, {-13.2343869456808, 517.908641346145}, {-13.1726447645511, 517.82374300373}});
    // std::cout << approx_equal(curve1[curve1.size()-2], curve1[curve1.size()-3]) << "\n";
    // Curve curve_1  = Curve({{-0.2320513169, -0.3384415771}, {-0.2061375252,-0.17794736420000001}});
    // Curve curve_2 = Curve({{-0.2320513169, -0.3384415771}, {-0.2061375252, -0.17794736420000001}});
    // Curve curve_1 = Curve({{0, 0}, {1, 0}});
    // Curve curve_2 = Curve({{0, 0}, {0, 1}});
    // std::cout << IntegralFrechet(curve_1, curve_2, ParamMetric::L1, 1, nullptr)
    // .compute_matching()
    // .cost;

// {
//     auto start_time = std::chrono::high_resolution_clock::now();
//     auto ret = simplification::imai_iri::simplify(curves[0], 5);
//     auto end_time = std::chrono::high_resolution_clock::now();
//     auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
//     std::cout << "Time: " << duration << "ms\n";
//     std::cout << "Cost: " << ret.first << "\n";
//     std::cout << "Curve:";
//     for (const auto& p: ret.second.get_points())
//         std::cout << " " << p;
//     std::cout << "\n";
//     std::cout << "Length: " << ret.second.size() << std::endl;
// }
// {
//     auto start_time = std::chrono::high_resolution_clock::now();
//     DTW distance(curves[0], curves[1]);
//     auto end_time = std::chrono::high_resolution_clock::now();
//     auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
//     std::cout << "Time: " << duration << "ms\n";
//     std::cout << "Cost: " << distance.cost() << "\n";
//     std::cout << "Matching:";
//     for (const auto& p: distance.matching())
//         std::cout << " (" << p.first << ", " << p.second << ")";
//     std::cout << std::endl;
// }
    return 0;
}
