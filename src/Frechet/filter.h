#ifndef FILTER_H
#define FILTER_H

#include "geom.h"
#include "Curve.h"
#include "Frechet/certificate.h"

class Filter {
private:
    Certificate cert;
    Curve const* curve1_pt;
    Curve const* curve2_pt;
    distance_t distance;

public:
    Filter(Curve const& curve1, Curve const& curve2, distance_t dist) {
        this->curve1_pt = &curve1;
        this->curve2_pt = &curve2;
        this->distance = dist;
        cert.setCurves(&curve1, &curve2);
        cert.setDistance(distance);
    }

    Certificate const& getCertificate() {
        return cert;
    }

    bool bichromaticFarthestDistance();
    bool greedy();
    bool adaptiveGreedy(PointID& pos1, PointID& pos2);
    bool adaptiveSimultaneousGreedy();
    bool negative(PointID pos1, PointID pos2);

    static bool isPointTooFarFromCurve(Point fixed, Curve const& curve,
        distance_t distance);
    static bool isFree(Point const& fixed, Curve const& var_curve,
        PointID start, PointID end, distance_t distance);
    static bool isFree(Curve const& curve1, PointID start1, PointID end1,
        Curve const& curve2, PointID start2, PointID end2, distance_t distance);
    static void increase(std::size_t& step);
    static void decrease(std::size_t& step);
};
#endif
