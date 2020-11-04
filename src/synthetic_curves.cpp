#include "synthetic_curves.h"

#include <cmath>
#include <fstream>
#include <string>

synth::DisplacementGenerator::DisplacementGenerator(double lower_lim,
    double upper_lim, std::mt19937_64::result_type seed):
        rng(seed), lower(lower_lim), upper(upper_lim) {}

Point synth::DisplacementGenerator::generate_disp_vec() {
    constexpr auto pi = 3.14159265358979323846;
    std::uniform_real_distribution<> uniform_dist_angle(0.0, 2.0 * pi);
    std::uniform_real_distribution<> uniform_dist_length(lower, upper);
    double theta = uniform_dist_angle(rng);
    double length = uniform_dist_length(rng);
    return {length * std::cos(theta), length * std::sin(theta)};
}

std::vector<Curve> synth::generate_curves(Curve const& curve,
        std::size_t count) {
    DisplacementGenerator disp_gen(0.0, 1.0);
    std::vector<Curve> curves;

    for (std::size_t i = 0; i < count; ++i) {
        Curve new_curve;
        auto translation = disp_gen.generate_disp_vec();
        for (auto const& p: curve.get_points()) {
            auto disp_vec = disp_gen.generate_disp_vec() / 20.0;
            Point new_point = p + translation + disp_vec;
            new_curve.push_back(new_point);
        }
        curves.emplace_back(new_curve);
    }
    return curves;
}

void synth::write_curves(Curve const& base) {
    auto synthetic_curves = generate_curves(base, 20);
    std::ofstream dataset("synthetic_curves/dataset.txt",
        std::fstream::out | std::fstream::trunc);

    for (std::size_t i = 0; i < synthetic_curves.size(); ++i) {
        std::string file_name = "synthetic_curves/curve_" + std::to_string(i)
            + ".txt";
        io::export_points(file_name, synthetic_curves[i].get_points());
        dataset << "curve_" + std::to_string(i) + ".txt" << "\n";
    }

    dataset.close();
    io::export_points("synthetic_curves/true_center.txt", base.get_points());
}
