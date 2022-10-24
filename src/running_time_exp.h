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
        
        for (int i = 0; i < 5; ++i) {

            auto start_time = std::chrono::high_resolution_clock::now();
            auto cdtw = CDTW<2, Norm::L1, Norm::L1>(curves[i], curves[i+1]);
            auto end_time = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
            total_time += duration;
            cdtw.find_max_pieces();
            piece_count += cdtw.max_pieces;
        }

        double p_av = piece_count / 10;
        double t_av = total_time / 10;

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

    void experiment() {
        auto curves = io::read_curves("data/characters/data");
        // auto curves = io::read_curves("data/Data_for_Mann_et_al_RSBL"
            // "/Bladon & Church route recapping/bladon heath/utm", 1);

        int seed = 10;
        std::mt19937_64 g(seed);

        std::shuffle(curves.begin(), curves.end(), g);

        auto trunc_instances = std::vector<Curves>();
        auto sample_instances = std::vector<Curves>();

        auto trunc_results = std::ofstream("trunc_results_char.txt");
        auto sample_results = std::ofstream("sample_results.txt");

        trunc_results << "i,number_of_pieces,time\n";
        sample_results << "i,number_of_pieces,time\n";

        for (int i = 10; i <= 100; i+=10) {

            std::cout << i << std::endl;
            
            trunc_instances.push_back(Curves());
            sample_instances.push_back(Curves());

            for (auto curve: curves) {
                trunc_instances.back().push_back(shorten_curve(curve, TRUNCATE, i));
                sample_instances.back().push_back(shorten_curve(curve, SAMPLE, i));
            }

            auto piece_time_trunc = get_average_piece_count_and_time(trunc_instances.back());
            auto piece_time_sample = get_average_piece_count_and_time(sample_instances.back());

            trunc_results << i << "," << piece_time_trunc.first << "," << piece_time_trunc.second << "\n";
            sample_results << i << "," << piece_time_sample.first << "," << piece_time_sample.second << "\n";
        }

        trunc_results.close();    
    }
}