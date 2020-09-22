#include "clustering/plot_clustering.h"

#include <fstream>
#include <stdexcept>

void clustering::plot_clustering(Clustering const& clustering,
        std::vector<Curve> const& curves, std::string const& filename) {
    std::ofstream output_file(filename);
    if (!output_file.is_open())
        throw std::runtime_error("Failed to open file " + filename);
    output_file << clustering.size() << "\n";

    for (auto const& cluster: clustering) {
        output_file << cluster.curve_ids.size() << "\n";
        for (auto const& curve_id: cluster.curve_ids) {
            auto const& curve = curves[curve_id];
            for (std::size_t i = 0; i < curve.size(); ++i) {
                output_file << curve[i].x << " " << curve[i].y;
                output_file << (i == curve.size() - 1 ? "\n" : " ");
            }
        }

        auto const& cc = cluster.center_curve;
        for (std::size_t i = 0; i < cc.size(); ++i) {
            output_file << cc[i].x << " " << cc[i].y;
            output_file << (i == cc.size() - 1 ? "\n" : " ");
        }
    }
}
