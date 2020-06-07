#include "clustering/center_algs.h"

//#include "curve_simplification.h"
#include "DTW/dtw.h"
#include "Frechet/frechet_light.h"
#include "Frechet/frechet_matching.h"
#include "geom.h"
#include "IntegralFrechet/IntegralFrechet.h"
#include <limits>

namespace
{

std::vector<distance_t> get_redistributed_distances(Curve curve, 
        double lower_threshold=1./10, double upper_threshold=1./5) {

    Points new_points = Points();

    //indices of new points
    std::vector<size_t> indices = std::vector<size_t>();
    
    int point_budget = 0;

    // First pass over curve
    // Remove redundant points
    for (std::size_t i = 0; i < curve.get_points().size(); ++i) {
        if (i > 0 && i < curve.get_points().size() - 1) {

            Point a = new_points.back();
            Point b = curve[i];
            Point c = curve[i+1];

            Point dir_1 = Point(b.x - a.x, b.y - a.y);
            Point dir_2 = Point(c.x - b.x, c.y - b.y);

            distance_t theta = acute_angle(dir_1, dir_2);
            // std::cout << theta << " >= " << lower_threshold << "\n";
            if (theta >= lower_threshold) {
                if (!approx_equal(new_points.back(), curve[i])) {
                    new_points.push_back(curve[i]);
                    indices.push_back(i);
                    continue;
                }
            } else {
                // std::cout << "skipped a point\n";
            }

            point_budget++;

        } else {
            if (i == 0 || !approx_equal(new_points.back(), curve[i])) {
                new_points.push_back(curve[i]);
                indices.push_back(i);
            }
        }
    }

    // Second pass over the curve, add extra
    // points to sections with high curvature
    Points final_points = Points();
    std::vector<distance_t> final_distances = std::vector<distance_t>();
    size_t last_index = 0;

    std::size_t i = 0;

    while (i < new_points.size()) {
        if (i > 0 && i < new_points.size() - 1) {
            
            Point a = final_points.back();
            Point b = new_points[i];
            Point c = new_points[i+1];

            Point dir_1 = Point(b.x - a.x, b.y - a.y);
            Point dir_2 = Point(c.x - b.x, c.y - c.y);

            distance_t theta = acute_angle(dir_1, dir_2);

            if (theta > upper_threshold && point_budget > 0) {
                std::cout << "new point...\n";
                distance_t dist_1 = a.dist(b);
                distance_t dist_2 = b.dist(c);

                assert(dist_1 > 0);
                assert(dist_2 > 0);

                if (dist_1 > dist_2) {
                    final_distances.push_back(final_distances.back() + 9* dist_1 / 10);
                    final_distances.push_back(final_distances.back() + dist_1 / 10);
                    final_distances.push_back(final_distances.back() + dist_2);
                } else {
                    final_distances.push_back(final_distances.back() + dist_1);
                    final_distances.push_back(final_distances.back() + dist_2 / 10);
                    final_distances.push_back(final_distances.back() + 9 * dist_2 / 10);
                }

                final_points.push_back(c);
                point_budget--;

                // Curve temp_curve = Curve({a, b, c});
                // distance_t total_length = temp_curve.curve_length();
                // distance_t third = total_length / 3;

                // CPoint cpoint_1 = temp_curve.get_cpoint_after(third);
                // CPoint cpoint_2 = temp_curve.get_cpoint_after(2 * third);

                // Point new_point_1 = temp_curve.interpolate_at(cpoint_1);
                // Point new_point_2 = temp_curve.interpolate_at(cpoint_2);

                // final_points.push_back(new_point_1);
                // final_points.push_back(new_point_2);
                // final_points.push_back(c);

                // distance_t dist_a = curve.curve_length(last_index);
                // distance_t new_dist_1 = curve.curve_length(
                //  curve.get_cpoint_after(
                //      curve.curve_length(last_index) + third
                //  )
                // );
                // distance_t new_dist_2 = curve.curve_length(
                //  curve.get_cpoint_after(
                //      curve.curve_length(last_index) + 2 * third
                //  )
                // );
                // distance_t dist_c = curve.curve_length(indices[i+1]);


                // std::cout << "extra points!!";
                // final_distances.push_back(new_dist_1);
                // final_distances.push_back(new_dist_2);
                // final_distances.push_back(dist_c);


                // last_index = indices[i+1];
                // point_budget--;
            } else {
                final_points.push_back(b);
                final_points.push_back(c);
                last_index = indices[i+1];

                final_distances.push_back(curve.curve_length(indices[i]));
                final_distances.push_back(curve.curve_length(indices[i+1]));
            }

            i += 2;

        } else {
            if (final_distances.size() > 0) {
                assert(!approx_equal(final_distances.back(), curve.curve_length(indices[i])));
            }
            final_points.push_back(new_points[i]);
            final_distances.push_back(curve.curve_length(indices[i]));
            i += 1;
        }
    }

    // for (auto dist: final_distances) {
    //  std::cout << dist << "\n";
    // }
    // std::cout << Points(curve.get_points()) << "\n";
    // std::cout << "length: " << curve.curve_length();
    return final_distances;
}

// Finds points on curve_2 alined with
// the points on another curve defined by distances.
// The param_space_path argument is a sequence of 
// points defining the matching between the two curves.

Points get_images_of_points(
    Points param_space_path, std::vector<distance_t> distances,  Curve curve_2
) {
    Points points = Points();

    std::vector<distance_t> x_coords = std::vector<distance_t>();
    for (auto p: param_space_path) {
        x_coords.push_back(p.x);
    }
    points.push_back(curve_2[0]);


    std::size_t last_index = 0;

    for (std::size_t i = 0; i < distances.size(); ++i) {

        distance_t distance = distances[i];
        
        if (distance == 0)
            continue;

        while (last_index <= x_coords.size() - 2 && x_coords[last_index + 1] <= distance) {
            ++last_index;
        }

        auto prev = last_index;
        auto next = last_index + 1;

        if (prev >= param_space_path.size()) {
            prev = param_space_path.size() - 2;
            next = prev + 1;
        }

        Point p = param_space_path[prev];
        Point q = param_space_path[next];

        Line line = Line::fromTwoPoints(p, q);
        distance_t distance_along_c2;
        if (line.isVertical()) {
            distance_along_c2 = p.y;
        } else {
            distance_along_c2 = line.getY(distances[i]);
            // std::cout << distance_along_c2 << ", " << curve_2.curve_length() << "\n";
            assert(distance_along_c2 >= 0);
        }

        if (distance_along_c2 > curve_2.curve_length()) {
            distance_along_c2 = curve_2.curve_length();
        }
        CPoint cpoint = curve_2.get_cpoint_after(distance_along_c2);
        Point new_point = curve_2.interpolate_at(cpoint);

        // if (i == distances.size() - 1) {
        //  if (!approx_equal(new_point, curve_2.get_points().back())) {
        //      std::cout << "uh oh...\n";
        //  }
        //  // assert(approx_equal(new_point, curve_2.get_points().back()));
        // }

        points.push_back(
            new_point
        );
    }

    return points;
}

std::vector<distance_t> get_uniform_distances_along_curve(Curve curve) {

    distance_t distance = curve.curve_length() / curve.size();

    return std::vector<distance_t>(curve.size() - 1, distance);
}

Points get_uniformly_spaced_points(Points param_space_path, Curve curve_1, Curve curve_2) {

    //sequence of points on curve_2
    Points points = Points();

    distance_t curve_1_length = curve_1.curve_length();
    distance_t interval_size = curve_1_length / (curve_1.get_points().size() - 1);

    std::vector<distance_t> x_coords = std::vector<distance_t>();
    for (auto p: param_space_path) {
        x_coords.push_back(p.x);
    }
    points.push_back(curve_2[0]);
    for (std::size_t i = 1; i < curve_1.get_points().size(); ++i) {
        auto it =  std::lower_bound(x_coords.begin(), x_coords.end(), interval_size * i);
        auto prev_index = static_cast<std::size_t>(std::distance(x_coords.begin(), it)) - 1;
        auto next_index = prev_index+1;
        
        if (prev_index >= param_space_path.size()) {
            prev_index = param_space_path.size() - 2;
            next_index = prev_index + 1;
        }

        Point p = param_space_path[prev_index];
        Point q = param_space_path[next_index];

        Line line = Line::fromTwoPoints(p, q);
        distance_t distance_along_c2;
        if (line.isVertical()) {
            distance_along_c2 = p.y;
        } else {
            distance_along_c2 = line.getY(interval_size * i);
            assert(distance_along_c2 > 0);
        }

        if (distance_along_c2 > curve_2.curve_length()) {
            distance_along_c2 = curve_2.curve_length();
        }
        CPoint cpoint = curve_2.get_cpoint_after(distance_along_c2);
        Point new_point = curve_2.interpolate_at(cpoint);
        
        points.push_back(
            new_point
        );
    }

    return points;
}


Point mean_of_points(Points points) {

    distance_t x_mean = 0;
    distance_t y_mean = 0;

    for (auto p: points) {
        x_mean += p.x;
        y_mean += p.y;
    }

    x_mean /= points.size();
    y_mean /= points.size();

    Point new_point = Point(x_mean, y_mean);

    return new_point;
}

bool calcKXCenters(clustering::Curves const& curves, Clustering& clustering,
        std::size_t l, clustering::C2CDist c2c_dist,
        std::function<distance_t(Curve const&, Curve const&)> const& dist) {
    bool found_new_center = false;

    // compute cluster costs in case they haven't been computed yet
    for (auto& cluster: clustering) {
        if (cluster.cost == std::numeric_limits<distance_t>::max()) {
            cluster.cost = calcC2CDist(curves, cluster.center_curve,
                cluster.curve_ids, c2c_dist, dist);
        }
    }

    for (auto& cluster: clustering) {
        for (CurveID curve_id1: cluster.curve_ids) {
            // auto simplified_curve = simplify(curves[curve_id1], l);
            auto simplified_curve = curves[curve_id1].naive_l_simplification(l);
            auto d = calcC2CDist(curves, simplified_curve, cluster.curve_ids,
                c2c_dist, dist);
            if (d < cluster.cost) {
                cluster.center_curve = simplified_curve;
                cluster.cost = d;
                found_new_center = true;
            }
        }
    }
    return found_new_center;
}

} // end anonymous namespace

std::string clustering::toString(CenterAlg center_alg) {
    switch(center_alg) {
    case CenterAlg::kMedian: return "kMedian";
    case CenterAlg::kMeans: return "kMeans";
    case CenterAlg::kCenter: return "kCenter";
    case CenterAlg::fCenter: return "fCenter";
    case CenterAlg::fMean: return "fMean";
    case CenterAlg::newCenterUpdate: return "newCenterUpdate";
    case CenterAlg::newCenterUpdate2: return "newCenterUpdate2";
    case CenterAlg::ensembleMethod1: return "ensembleMethod1";
    }
    ERROR("Unknown center_alg.");
}

distance_t clustering::calcC2CDist(Curves const& curves,
        Curve const& center_curve, CurveIDs const& curve_ids, C2CDist c2c_dist,
        std::function<distance_t(Curve const&, Curve const&)> const& dist) {
    distance_t d = 0.0;
    for (auto const& curve_id: curve_ids) {
        auto curve_dist = dist(center_curve, curves[curve_id]);

        switch (c2c_dist) {
        case C2CDist::Median:
            d += curve_dist;
            break;
        case C2CDist::Mean:
            d += curve_dist * curve_dist;
            break;
        case C2CDist::Max:
            d = std::max(d, curve_dist);
            break;
        }
    }
    return d;
}

bool clustering::computerCenters(Curves const& curves, Clustering& clustering,
        std::size_t l, CenterAlg center_alg,
        std::function<distance_t(Curve const&, Curve const&)> const& dist) {
    switch (center_alg) {
    case CenterAlg::kMedian:
        return calcKMedianCenters(curves, clustering, l, dist);
    case CenterAlg::kMeans:
        return calcKMeansCenters(curves, clustering, l, dist);
    case CenterAlg::kCenter:
        return calcKCenterCenters(curves, clustering, l, dist);
    case CenterAlg::fCenter:
        return calcFSACenters(curves, clustering, l, dist, C2CDist::Max,
            CenterCurveUpdateMethod::frechetCentering);
    case CenterAlg::fMean:
        return calcFSACenters(curves, clustering, l, dist, C2CDist::Median,
            CenterCurveUpdateMethod::frechetMean);
    case CenterAlg::dtwMean:
        return calcFSACenters(curves, clustering, l, dist, C2CDist::Median,
            CenterCurveUpdateMethod::dtwMean);
    case CenterAlg::avFCenter:
        return calcFSACenters(curves, clustering, l, dist, C2CDist::Median,
            CenterCurveUpdateMethod::avFCenter);
    case CenterAlg::newCenterUpdate:
        return calcFSACenters(curves, clustering, l, dist, C2CDist::Median,
            CenterCurveUpdateMethod::newCenterUpdate);
    case CenterAlg::newCenterUpdate2:
        return calcFSACenters(curves, clustering, l, dist, C2CDist::Median,
            CenterCurveUpdateMethod::newCenterUpdate2);
    case CenterAlg::naiveCenterUpdate:
        return naiveCenterUpdate(curves, clustering, dist, C2CDist::Median);
    case CenterAlg::ensembleMethod1:
        return ensembleMethod1(curves, clustering, dist, C2CDist::Median);
    }
    ERROR("No matching center_alg enum passed.");
}

bool clustering::calcKMedianCenters(Curves const& curves,
        Clustering& clustering, std::size_t l,
        std::function<distance_t(Curve const&, Curve const&)> const& dist) {
    return calcKXCenters(curves, clustering, l, C2CDist::Median, dist);
}

bool clustering::calcKMeansCenters(Curves const& curves,
        Clustering& clustering, std::size_t l,
        std::function<distance_t(Curve const&, Curve const&)> const& dist) {
    return calcKXCenters(curves, clustering, l, C2CDist::Mean, dist);
}

bool clustering::calcKCenterCenters(Curves const& curves,
        Clustering& clustering, std::size_t l,
        std::function<distance_t(Curve const&, Curve const&)> const& dist) {
    return calcKXCenters(curves, clustering, l, C2CDist::Max, dist);
}

// bool frechetCentering (Curves const& curves, Clustering& clustering, int l, C2CDist c2c_dist, distance_t(*dist_func)(Curve, Curve)) {
//     return calcFSACenters(curves, clustering, l, dist_func, c2c_dist, CenterCurveUpdateMethod::frechetCentering);
// }

// bool frechetMean(Curves const& curves, Clustering& clustering, int l, C2CDist c2c_dist, distance_t(*dist_func)(Curve, Curve)) {
//     return calcFSACenters(curves, clustering, l, dist_func, c2c_dist, CenterCurveUpdateMethod::frechetMean);
// }

// TODO: There is some unnecessary pushing around of data here. Fix that to increase performance.
bool clustering::calcFSACenters(Curves const& curves, Clustering& clustering,
        std::size_t /* l */,
        std::function<distance_t(Curve const&, Curve const&)> const& dist,
        C2CDist c2c_dist, CenterCurveUpdateMethod method) {
    bool found_new_center = false;

    // compute cluster costs in case they haven't been computed yet
    for (auto& cluster: clustering) {
        if (cluster.cost == std::numeric_limits<distance_t>::max()) {
            cluster.cost = calcC2CDist(curves, cluster.center_curve,
                cluster.curve_ids, c2c_dist, dist);
        }
    }

    for (auto& cluster: clustering) {
        std::vector<Points> matchings;
        auto const& center_curve = cluster.center_curve;
        Curve new_center_curve;

        for (auto curve_id: cluster.curve_ids) {
            auto const& curve = curves[curve_id];
            // auto matching = calcMatching(cluster.center_curve, curve);
            // auto matching = IntegralFrechet(cluster.center_curve, curve, ParamMetric::LInfinity_NoShortcuts, 1, nullptr)
            // .compute_matching()
            // .matching;

            // auto matching = method == frechetCentering ? 
            // calcMatching(cluster.center_curve, curve) :
            // IntegralFrechet(cluster.center_curve, curve, ParamMetric::LInfinity_NoShortcuts, 1, nullptr)
            // .compute_matching()
            // .matching;

            // Points matching_points = matching_to_points(matching, cluster.center_curve, curve);
            // matchings.push_back(std::move(matching_points));

            switch (method) {
                case CenterCurveUpdateMethod::frechetCentering:
                    matchings.push_back(calcMatching(cluster.center_curve, curve));
                    break;
                case CenterCurveUpdateMethod::avFCenter:
                    matchings.push_back(
                        get_images_of_points(
                            IntegralFrechet(cluster.center_curve, curve, ParamMetric::LInfinity_NoShortcuts, 1, nullptr)
                            .compute_matching()
                            .matching,
                            cluster.center_curve.get_prefix_length_vector(),
                            curve
                        )
                    );
                    break;
                case CenterCurveUpdateMethod::frechetMean:
                    matchings.push_back(
                        get_images_of_points(
                            IntegralFrechet(cluster.center_curve, curve, ParamMetric::LInfinity_NoShortcuts, 1, nullptr)
                            .compute_matching()
                            .matching,
                            cluster.center_curve.get_prefix_length_vector(),
                            curve
                        )
                    );
                    break;
                case CenterCurveUpdateMethod::dtwMean: {
                    auto dtw_matching = DTW(cluster.center_curve, curve).matching();
                    Points matching = Points();
                    std::size_t index = 0;
                    for (std::size_t i = 0; i < cluster.center_curve.get_points().size(); ++i) {
                        while (i != dtw_matching[index].first) {
                            index++;
                        }
                        matching.push_back(curve[index]);
                    }
                    break;
                }
                case CenterCurveUpdateMethod::newCenterUpdate:
                    matchings.push_back(
                        get_images_of_points(
                            IntegralFrechet(cluster.center_curve, curve, ParamMetric::LInfinity_NoShortcuts, 1, nullptr)
                            .compute_matching()
                            .matching,
                            get_uniform_distances_along_curve(cluster.center_curve),
                            curve
                        )
                    );
                    break;
                case CenterCurveUpdateMethod::newCenterUpdate2:
                    matchings.push_back(
                        get_images_of_points(
                            IntegralFrechet(cluster.center_curve, curve, ParamMetric::LInfinity_NoShortcuts, 1, nullptr)
                            .compute_matching()
                            .matching,
                            get_redistributed_distances(cluster.center_curve),
                            curve
                        )
                    );
                    break;
            }
        }


        for (PointID point_id = 0; point_id < center_curve.size(); ++point_id) {
            Points matching_points;
            for (auto const& matching: matchings) {
                matching_points.push_back(matching[point_id]);
            }

            Point new_point;

            switch (method) {
                case CenterCurveUpdateMethod::frechetCentering:
                    new_point = calcMinEnclosingCircle(matching_points).center;
                    break;
                case CenterCurveUpdateMethod::avFCenter:
                    new_point = calcMinEnclosingCircle(matching_points).center;
                    break;
                case CenterCurveUpdateMethod::frechetMean:
                    new_point = mean_of_points(matching_points);
                    break;
                case CenterCurveUpdateMethod::dtwMean:
                    new_point = mean_of_points(matching_points);
                    break;
                case CenterCurveUpdateMethod::newCenterUpdate:
                    new_point = mean_of_points(matching_points);
                    break;
                case CenterCurveUpdateMethod::newCenterUpdate2:
                    new_point = mean_of_points(matching_points);
                    break;
            }

            if (new_center_curve.get_points().empty() || !approx_equal(new_point, new_center_curve.get_points().back())) {
                new_center_curve.push_back(new_point);
            }
        }

        if (center_curve != new_center_curve) {
            auto new_dist = calcC2CDist(curves, new_center_curve, cluster.curve_ids, c2c_dist, dist);
            if (new_dist < cluster.cost) {
                cluster.center_curve = std::move(new_center_curve);
                cluster.cost = new_dist;
                found_new_center = true;
            }
        } 
    }
    
    if (found_new_center) {
        std::cout << "found new center\n";
    }
    else {
        std::cout << "no new center... :( \n";
    }
    return found_new_center;
}

bool clustering::naiveCenterUpdate(Curves const& curves, Clustering& clustering,
        std::function<distance_t(Curve const&, Curve const&)> const& dist,
        C2CDist c2c_dist) {
    bool found_new_center = false;

    for (auto& cluster: clustering) {
        if (cluster.cost == std::numeric_limits<distance_t>::max()) {
            cluster.cost = calcC2CDist(curves, cluster.center_curve, cluster.curve_ids, c2c_dist, dist);
        }
    }

    for (auto& cluster: clustering) {
        auto& center_curve = cluster.center_curve;
        auto l = center_curve.get_points().size();
        std::vector<Points> matchings = std::vector<Points>();
        for (std::size_t i = 0; i < l; ++i) {
            matchings.push_back(Points());
        }


        for (auto& curve_id: cluster.curve_ids) {
            Curve curve = curves[curve_id];
            distance_t d = curve.curve_length() / (l - 1);

            for (std::size_t i = 0; i < l; ++i) {
                matchings[i].push_back(curve.interpolate_at(curve.get_cpoint_after(d * i)));
            }
        }
        std::cout << "finished first loop\n";
        Points new_points = Points();

        for (std::size_t i = 0; i < matchings.size(); ++i) {
            if (matchings[i].size() > 0) {
                Point point = mean_of_points(matchings[i]);
                // if (new_points.size() > 0 && !approx_equal(new_points.back(), point))
                new_points.push_back(point);
            }
        }

        Curve new_center_curve = Curve(new_points);

        if (center_curve != new_center_curve) {
            auto new_dist = calcC2CDist(curves, new_center_curve, cluster.curve_ids, c2c_dist, dist);
            if (new_dist < cluster.cost) {
                cluster.center_curve = std::move(new_center_curve);
                cluster.cost = new_dist;
                found_new_center = true;
            }
        } 

    }

    if (found_new_center) {
        std::cout << "found new center\n";
    }
    else {
        std::cout << "no new center... :( \n";
    }
    return found_new_center;

}

bool clustering::ensembleMethod1(Curves const& curves, Clustering& clustering,
        std::function<distance_t(Curve const&, Curve const&)> const& dist,
        C2CDist c2c_dist) {

    if (naiveCenterUpdate(curves, clustering, dist, c2c_dist)) {
        std::cout << "naive method helped...\n";
    }

    bool new_centers = calcFSACenters(curves, clustering, 0, dist, c2c_dist, CenterCurveUpdateMethod::frechetMean);
    return new_centers;
}
