#pragma once

#include "distance_functions.h"
#include "io.cpp"
#include <fstream>
#include <math.h>

class DisplacementGenerator {
    
    std::mt19937_64 rng;
    double lower_lim;
    double upper_lim;
    
    public:
        DisplacementGenerator(int seed, double lower_lim, double upper_lim) : 
        rng(std::mt19937_64(seed)), lower_lim(lower_lim), upper_lim(upper_lim) {};

        Point generate_disp_vec() {
            
            std::uniform_real_distribution<double> uniform_dist_angle(0, 2*M_PI);
            std::uniform_real_distribution<double> uniform_dist_length(lower_lim, upper_lim);
            double theta = uniform_dist_angle(rng);
            double length = uniform_dist_length(rng);
            
            return {length * std::cos(theta), length * std::sin(theta)};
        };
};

Curves generate_curves(Curve curve, int count) {

    DisplacementGenerator disp_gen = DisplacementGenerator(1234, 0, 1);
    
    Curves curves = Curves();

    for (int i = 0; i < count; ++i) {
        Curve new_curve = Curve();
        Point translation = disp_gen.generate_disp_vec();
        for (auto& p: curve.get_points()) {
            Point disp_vec = disp_gen.generate_disp_vec() / 20;
            Point new_point = {p.x + translation.x + disp_vec.x, p.y + translation.y + disp_vec.y};

            new_curve.push_back(new_point);
        }

        // std::cout << new_curve.get_points() << "\n";
        curves.push_back(new_curve);
    }


    return curves;
}

void write_curves() {
    Curves curves = io::read_curves("data/characters/data");
    Curve curve = curves[0];
    
    Curves synthetic_curves = generate_curves(curve, 20);
    std::fstream dataset;
    dataset.open("synthetic_curves/dataset.txt", std::fstream::out | std::fstream::trunc);

    for (int i = 0; i < synthetic_curves.size(); ++i) {
        Curve synth_curve = synthetic_curves[i];
        std::string file_name = "synthetic_curves/curve_" + std::to_string(i) + ".txt";
        io::export_points(file_name, synth_curve.get_points());
        dataset << "curve_" + std::to_string(i) + ".txt" << "\n";
    }

    dataset.close();
    io::export_points("synthetic_curves/true_center.txt", curve.get_points());
}

