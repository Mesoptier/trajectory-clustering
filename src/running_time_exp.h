#include "cdtw/2d-l1-l1.h"
#include "cdtw/cdtw.h"
#include <algorithm>
#include <random>


namespace running_time {

    enum SCALING_METHOD {
        TRUNCATE, SAMPLE
    };

    std::pair<double, double> get_average_piece_count_and_time(Curves curves) {
        double piece_count = 0;
        double total_time = 0;

        int number_of_smaples = 5;
        
        for (int i = 0; i < number_of_smaples; ++i) {
            auto start_time = std::chrono::high_resolution_clock::now();
            auto cdtw = CDTW<2, Norm::L1, Norm::L1>(curves[i], curves[i+1]);
            auto end_time = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
            total_time += duration;
            cdtw.find_max_pieces();
            piece_count += cdtw.max_pieces;
        }

        double p_av = piece_count / number_of_smaples;
        double t_av = total_time / number_of_smaples;

        std::pair<double, double> result(p_av, t_av);

        return result;
    }


    Curve shorten_curve(const Curve& curve, SCALING_METHOD method, int size) {
        switch(method) {
            case TRUNCATE: {
                return curve.slice(0, size);
                break;
            }
            case SAMPLE: {
                return curve.naive_l_simplification(size);
                break;
            }
        }
    }

    void experiment(std::string dataset) {
        Curves curves;
        if (dataset == "pigeons") {
            curves = io::read_curves("data/Data_for_Mann_et_al_RSBL"
            "/Bladon & Church route recapping/bladon heath/utm", 1);
        } else if (dataset == "characters") {
            curves = io::read_curves("data/characters/data");
        }

         int seed = 10;
        std::mt19937_64 g(seed);

        std::shuffle(curves.begin(), curves.end(), g);

        auto instances = std::vector<Curves>();

        auto results = std::ofstream("2d_paper_results/running_time_results_" + dataset + ".csv");

        results << "i,number_of_pieces,time\n";

        for (int i = 10; i <= 100; i+=10) {

            std::cout << i << std::endl;
            
            instances.push_back(Curves());

            for (auto curve: curves) {
                instances.back().push_back(shorten_curve(curve, TRUNCATE, i));
            }

            auto piece_time_trunc = get_average_piece_count_and_time(instances.back());

            results << i << "," << piece_time_trunc.first << "," << piece_time_trunc.second << "\n";
        }

        results.close();
    }

    void run_experiments() {
        experiment("pigeons");
        experiment("characters");
    }
}