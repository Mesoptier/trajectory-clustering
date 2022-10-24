#include "distance_functions.h"

#include "DTW/dtw.h"
#include "Frechet/frechet_light.h"
#include "IntegralFrechet/IntegralFrechet.h"
#include "cdtw/cdtw.h"
#include "cdtw/2d-l1-l1.h"

namespace {
    /**
     * \brief Reparametrize a matching over simplified curves to a matching over
     * the original curves. (See also in main.cpp.)
     *
     * \param matching - Simplified matching, will be reparametrized in-place.
     * \param c1 - Original curve for x-axis
     * \param c1_s - Simplified curve for x-axis
     * \param c2 - Original curve for y-axis
     * \param c2_s - Simplified curve for y-axis
     */
    void reparametrize_matching(Points& matching,
            Curve const& c1, SimplifiedCurve const& c1_s,
            Curve const& c2, SimplifiedCurve const& c2_s) {
        PointID prev_id1 = 0;
        PointID prev_id2 = 0;

        for (auto& p : matching) {
            // For p.x/p.y:
            // - Find PointID's of corresponding edge in (c1_s/c2_s).curve
            // - Find corresponding PointID's in c1/c2 from
            //   (c1_s/c2_s).original_points
            // - Update p.x/p.y by linear interpolation

            auto const cp1 = c1_s.curve.get_cpoint_after(p.x, prev_id1);
            prev_id1 = cp1.getPoint();
            p.x = c1.curve_length(c1_s.original_points[cp1.getPoint()])
                * (1 - cp1.getFraction());
            if (cp1.getFraction() != 0) {
                p.x += c1.curve_length(c1_s.original_points[cp1.getPoint() + 1])
                    * cp1.getFraction();
            }

            auto const cp2 = c2_s.curve.get_cpoint_after(p.y, prev_id2);
            prev_id2 = cp2.getPoint();
            p.y = c2.curve_length(c2_s.original_points[cp2.getPoint()])
                * (1 - cp2.getFraction());
            if (cp2.getFraction() != 0) {
                p.y += c2.curve_length(c2_s.original_points[cp2.getPoint() + 1])
                    * cp2.getFraction();
            }
        }
    }
}

distance_t df::dtw(Curve const& curve_1, Curve const& curve_2) {
    return DTW(curve_1, curve_2).cost();
}

void write_heur_warping_path(IntegralFrechet::MatchingResult matching) {
    std::ofstream file("warping_path_h.txt");

    auto m = matching.matching;

    for (int i = 1; i < m.size(); ++i) {
        file << m[i-1].x << " " << m[i-1].y << " "
        << m[i].x << " " << m[i].y << "\n";
    }

    file.close();
}

distance_t df::integral_frechet(Curve const& curve_1, Curve const& curve_2) {
    auto const res = std::max(static_cast<std::size_t>(
        (curve_1.curve_length() + curve_2.curve_length())
        / (curve_1.size() + curve_2.size()) / 5), 1UL);

    return IntegralFrechet(curve_1, curve_2, ParamMetric::L1, res, nullptr)
        .compute_matching().cost;
}

distance_t df::integral_frechet_fast(Curve const& curve_1,
        Curve const& curve_2) {
    auto const res = std::max(static_cast<std::size_t>(
        (curve_1.curve_length() + curve_2.curve_length())
        / (curve_1.size() + curve_2.size()) / 5), 1UL);
    auto const curve1_simpl = curve_1.simplify();
    auto const curve2_simpl = curve_2.simplify();
    if (curve1_simpl.curve == curve2_simpl.curve)
        return integral_frechet(curve_1, curve_2);
    auto matching_simpl = IntegralFrechet(curve1_simpl.curve,
        curve2_simpl.curve, ParamMetric::L1, 10 * res, nullptr, ImageMetric::L2)
        .compute_matching().matching;
    reparametrize_matching(matching_simpl, curve_1, curve1_simpl,
        curve_2, curve2_simpl);

    MatchingBand const band(curve_1, curve_2, matching_simpl, 5 * res);
    return IntegralFrechet(curve_1, curve_2, ParamMetric::L1, res, &band)
        .compute_matching().cost;
}

distance_t df::frechet(Curve const& curve_1, Curve const& curve_2) {
    FrechetLight fl;
    return fl.calcDistance(curve_1, curve_2);
}

distance_t df::average_frechet(Curve const& curve_1, Curve const& curve_2) {
    return integral_frechet(curve_1, curve_2) /
        (curve_1.curve_length() + curve_2.curve_length());
}

distance_t df::average_frechet_fast(Curve const& curve_1,
        Curve const& curve_2) {
    return integral_frechet_fast(curve_1, curve_2) /
        (curve_1.curve_length() + curve_2.curve_length());
}

distance_t df::cdtw_2d_l1_l1(Curve const& curve_1,
    Curve const& curve_2) {

        using CDTW = CDTW<2, Norm::L1, Norm::L1>;
        
        auto cdtw = CDTW(curve_1, curve_2);
        return cdtw.cost();
}

distance_t df::heur_cdtw_2d_l1_l1(Curve const& curve_1,
    Curve const& curve_2) {
        auto const res = std::max(static_cast<std::size_t>(
        (curve_1.curve_length() + curve_2.curve_length())
        / (curve_1.size() + curve_2.size()) / 5), 1UL);
    write_heur_warping_path(IntegralFrechet(curve_1, curve_2, ParamMetric::L1, res, nullptr, ImageMetric::L1)
        .compute_matching());
    return IntegralFrechet(curve_1, curve_2, ParamMetric::L1, res, nullptr, ImageMetric::L1)
        .compute_matching().cost;
}

bool df::dtw_lt(Curve const& curve_1, Curve const& curve_2, distance_t delta) {
    return dtw(curve_1, curve_2) <= delta;
}

bool df::integral_frechet_lt(Curve const& curve_1, Curve const& curve_2,
        distance_t delta) {
    return integral_frechet(curve_1, curve_2) <= delta;
}

bool df::integral_frechet_fast_lt(Curve const& curve_1, Curve const& curve_2,
        distance_t delta) {
    return integral_frechet_fast(curve_1, curve_2) <= delta;
}

bool df::frechet_lt(Curve const& curve_1, Curve const& curve_2,
        distance_t delta) {
    FrechetLight fl;
    return fl.lessThanWithFilters(delta, curve_1, curve_2);
}

bool df::average_frechet_lt(Curve const& curve_1, Curve const& curve_2,
        distance_t delta) {
    return average_frechet(curve_1, curve_2) <= delta;
}

bool df::average_frechet_fast_lt(Curve const& curve_1, Curve const& curve_2,
        distance_t delta) {
    return average_frechet_fast(curve_1, curve_2) <= delta;
}
