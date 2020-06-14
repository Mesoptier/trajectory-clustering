#pragma once

#include "../basic_types.h"
#include "../defs.h"
#include <fstream>
#include <iostream>


void plot_clustering(Clustering const& clustering, Curves const& curves, std::string filename);

void plot_clustering(Clustering const& clustering, Curves const& curves, std::string filename) {

    std::fstream output_file;
    output_file.open(filename, std::fstream::out | std::fstream::trunc);

    output_file << clustering.size() << "\n";

    for (auto& cluster: clustering) {
        output_file << cluster.curve_ids.size() << "\n";
        for (CurveID curve_id: cluster.curve_ids) {
            Curve curve = curves[curve_id];
            for (Point point: curve.get_points()) {
                output_file << point.x << " ";
                if (point == curve.back())
                    output_file << point.y << "\n";
                else
                    output_file << point.y << " ";
            }
        }

        for (Point point: cluster.center_curve.get_points()) {
            output_file << point.x << " ";
            if (point == cluster.center_curve.back())
                output_file << point.y << "\n";
            else
                output_file << point.y << " ";
        }
    }
}