#ifndef CERTIFICATE_H
#define CERTIFICATE_H

#include "geom.h"
#include "Curve.h"

class Certificate {
public:
    bool isYes() const { 
        assert(isValid());
        return lessThan;
    }

    bool isValid() const {
        return valid;
    }

    bool check() const;
    CPositions const& getTraversal() const {
        return traversal;
    }

    void addPoint(CPosition const& pos) {
        traversal.push_back(pos);
    }

    void setAnswer(bool answer) {
        lessThan = answer;
    }

    void setCurves(Curve const* curve1, Curve const* curve2) {
        curve_pair[0] = curve1;
        curve_pair[1] = curve2;
    }

    void setDistance(distance_t distance) {
        dist = distance;
        dist_sqr = distance * distance;
    }

    void validate() {
        valid = true;
    }

    void reset() {
        valid = false;
        traversal.clear();
    }

    void clear() {
        traversal.clear();
    }

    void dump_certificate() const;

private:
    // if YES certificate: (traversal1[0], traversal2[0]), ...,
    // (traversal1[T], traversal2[T]) should be feasible traversal of
    // curve1 and curve2 (interpolate between discrete points).
    // if NO certificate: certificate should be as described in
    // certificate-outline.pdf

    CPositions traversal;
    std::array<Curve const*, 2> curve_pair;
    distance_t dist, dist_sqr;

    bool lessThan;
    bool valid = false;

    bool feasible(CPosition const& pt) const;
    bool feasible(CPoint const& pt1, CPoint const& pt2) const;
    bool nonEmpty(CurveID fixed_curve, CPoint const& fixed_point,
        CPoint const& start_pt, CPoint const& end_point) const;
};
#endif
