#pragma once

#include "../Curve.h"
// #include "../greedy_simplification.h"
#include <limits>

using Curves = std::vector<Curve>;


class Clustering {

    private:
        int k;
        int l;
        Curves curves;
        Curves centers;
        std::vector<std::vector<int>> clusters;


    public:
        Clustering(Curves curves, int k, int l) : curves(curves), k(k), l(l) {
            clusters = std::vector<std::vector<int>>();
            centers = Curves();

            for (int i = 0; i < k; ++i) {
                //std::cout << i << "\n";
                centers.push_back(
                    greedy_simplification(curves[i], 0.25)
                );

                clusters.push_back(
                    std::vector<int>()
                );
            }

            std::cout << clusters.size() << "\n";
        };

        void assign_centers() {
            for (int i = 0; i < clusters.size(); ++i) {
                clusters[i] = std::vector<int>();
            }

            for (int j = 0; j < curves.size(); ++j) {

                int min_index = -1;
                distance_t min_distance = std::numeric_limits<double>::max();

                for (int i = 0; i < k; ++i) {
                    
                    distance_t dist = IntegralFrechet(curves[j], centers[i], ParamMetric::LInfinity_NoShortcuts, 50, nullptr)
                    .compute_matching()
                    .cost;
                    
                    if (dist < min_distance) {
                        min_distance = dist;
                        min_index = i;
                    }

                }

                clusters[min_index].push_back(j);
            }
        }

        std::vector<std::vector<int>> get_clusters() {
            return clusters;
        }


        distance_t cost() {

            distance_t total_cost = 0;

            for (int i = 0; i < k; ++i) {
                for (int curve_id: clusters[i]) {
                    total_cost += IntegralFrechet(curves[curve_id], centers[i], ParamMetric::LInfinity_NoShortcuts, 1, nullptr)
                    .compute_matching()
                    .cost;
                }
            }

            return total_cost;
        }
};