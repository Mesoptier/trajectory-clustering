#pragma once

#include "../Curve.h"

class MatchingBand
{
    std::vector<CPoint> lower_x;
    std::vector<CPoint> upper_x;
    std::vector<CPoint> lower_y;
    std::vector<CPoint> upper_y;

public:
    MatchingBand(const Curve& curve_x, const Curve& curve_y, const Points& matching, distance_t radius);

    CPoint get_min_h(PointID) const;
    CPoint get_max_h(PointID) const;
    CPoint get_min_v(PointID) const;
    CPoint get_max_v(PointID) const;

};