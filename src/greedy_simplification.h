#pragma once
#include "Curve.h"
#include "IntegralFrechet/IntegralFrechet.h"





// Returns a simplification of the input curve following
// something resembling the greedy algorithm of Agarwal et al.
// Currently brute force so will be very slow...

Curve greedy_simplification(Curve const& curve, distance_t budget) {
    
    Curve simplification(curve.name());
    Points points = curve.get_points();
    std::vector<int> indices = {0};
    int last_index = 0;
    int size = points.size();

    while (last_index < size-1) {
        for (int i = last_index+1; i < size; ++i) {
            if (i == size - 1) {
                indices.push_back(i);
                last_index = size-1;
                continue;

            } else {

                Curve edge("edge", {points[last_index], points[i]});
                Curve next_edge("next_edge", {points[last_index], points[i+1]});
                Curve subcurve("subcurve", Points(points.begin() + last_index, points.begin() + i + 1));
                Curve next_subcurve("next_subcurve", Points(points.begin() + last_index, points.begin() + i + 2));
                

                distance_t edge_dist;
                if (subcurve.get_points().size() == edge.get_points().size()) {
                    edge_dist = 0;
                } else {
                    edge_dist = IntegralFrechet(edge, subcurve,  ParamMetric::LInfinity_NoShortcuts, 1, nullptr).compute_matching().cost;
                }

                distance_t next_edge_dist = IntegralFrechet(next_edge, next_subcurve, ParamMetric::LInfinity_NoShortcuts, 1, nullptr).compute_matching().cost; 
            
                if (edge_dist <= budget && next_edge_dist > budget) {
                    indices.push_back(i);
                    last_index = i;
                    continue;
                }

            }
        }
    }

    for (auto i: indices) {
        simplification.push_back(points[i]);
    }

    //distance_t final_distance = IntegralFrechet(simplification, curve, ParamMetric::L1, 100, nullptr).compute_matching().cost;
    //std::cout << final_distance / (curve.curve_length() + simplification.curve_length()) << " - final average distance" << "\n";

    return simplification;
}

// Curve greedy_l_simplification(Curve const& curve, int l) {

//     distance_t epsilon = 1e-8;
//     distance_t min = 0;
//     distance_t max = 100000000000;

//     Curve simplified_curve = curve;
// 	while (max-min > epsilon) {
// 		auto split = (max + min)/2.;
// 		simplified_curve = greedy_simplification(curve, split);
// 		if ((int)simplified_curve.size() <= l) {
// 			max = split;
// 		}
// 		else {
// 			min = split;
// 		}

//     return simplified_curve;
// }