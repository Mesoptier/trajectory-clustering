#ifndef CENTER_ALGS_H
#define CENTER_ALGS_H

#include <functional>
#include <string>
#include <vector>

#include "basic_types.h"

namespace clustering {
    using Curves = std::vector<Curve>;

    enum class CenterAlg {
        kMedian,
        kMeans,
        kCenter,
        fCenter,
        fMean,
        dtwMean,
        avFCenter,
        newCenterUpdate,
        newCenterUpdate2,
        naiveCenterUpdate,
        ensembleMethod1,
        dba,
        cdba,
        wedge
    };

    enum class CenterCurveUpdateMethod {
        frechetCentering,
        frechetMean,
        dtwMean,
        avFCenter,
        newCenterUpdate,
        newCenterUpdate2,
        naiveCenterUpdate,
        ensembleMethod1
    };

    std::string toString(CenterAlg center_alg);

    enum class C2CDist {
        Median,
        Mean,
        Max,
    };

    distance_t calcC2CDist(Curves const& curves, Curve const& center_curve,
        CurveIDs const& curve_ids, C2CDist c2c_dist,
        std::function<distance_t(Curve const&, Curve const&)> const& dist);

    bool computerCenters(Curves const& curves, Clustering& clustering,
        std::size_t l, CenterAlg center_alg,
        std::function<distance_t(Curve const&, Curve const&)> const& dist);

    bool calcKMedianCenters(Curves const& curves, Clustering& clustering,
        std::size_t l,
        std::function<distance_t(Curve const&, Curve const&)> const& dist);

    bool calcKMeansCenters(Curves const& curves, Clustering& clustering,
        std::size_t l,
        std::function<distance_t(Curve const&, Curve const&)> const& dist);

    bool calcKCenterCenters(Curves const& curves, Clustering& clustering,
        std::size_t l,
        std::function<distance_t(Curve const&, Curve const&)> const& dist);

    bool calcFSACenters(Curves const& curves, Clustering& clustering,
        std::size_t l,
        std::function<distance_t(Curve const&, Curve const&)> const& dist,
        C2CDist cluster_dist = C2CDist::Max,
        CenterCurveUpdateMethod method = CenterCurveUpdateMethod::frechetMean);

    bool naiveCenterUpdate(Curves const& curves, Clustering& clustering,
        std::function<distance_t(Curve const&, Curve const&)> const& dist,
        C2CDist c2c_dist);

    bool dba(Curves const& curves, Clustering& clustering, distance_t(*dist_func)(Curve, Curve), C2CDist c2c_dist);
    bool cdba(Curves const& curves, Clustering& clustering, distance_t(*dist_func)(Curve, Curve), C2CDist c2c_dist);

    // First computes a new center via the naiveCenterUpdateMethod
    // the frechetMean update method is then applied
    bool ensembleMethod1(Curves const& curves, Clustering& clustering,
        std::function<distance_t(Curve const&, Curve const&)> const& dist,
        C2CDist c2c_dist);
    // Points matching_of_vertices(Curve curve_1, Curve curve_2);
}
#endif
