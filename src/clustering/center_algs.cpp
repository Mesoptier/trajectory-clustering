#include "clustering/center_algs.h"

std::string clustering::toString(CenterAlg center_alg) {
    switch(center_alg) {
    case CenterAlg::fsa: return "FSA";
    case CenterAlg::dba: return "DBA";
    case CenterAlg::cdba: return "CDBA";
    case CenterAlg::wedge: return "Wedge method";
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
