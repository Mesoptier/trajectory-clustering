#pragma once

#include "../Curve.h"

class MatchingBand
{
    std::vector<CPoint> lower_y_at_x;
    std::vector<CPoint> upper_y_at_x;
    std::vector<CPoint> lower_x_at_y;
    std::vector<CPoint> upper_x_at_y;

public:
    MatchingBand(const Curve& curve_x, const Curve& curve_y, const Points& matching, distance_t radius);

    CPoint get_lower_y_at_x(PointID x) const { return lower_y_at_x.at(x); }
    CPoint get_upper_y_at_x(PointID x) const { return upper_y_at_x.at(x); }
    CPoint get_lower_x_at_y(PointID y) const { return lower_x_at_y.at(y); }
    CPoint get_upper_x_at_y(PointID y) const { return upper_x_at_y.at(y); }

    bool contains_point(const CPoint& x, const CPoint& y) const {
        if (x.getFraction() == 0) {
            return y >= get_lower_y_at_x(x.getPoint()) && y <= get_upper_y_at_x(x.getPoint());
        }
        if (y.getFraction() == 0) {
            return x >= get_lower_x_at_y(y.getPoint()) && x <= get_upper_x_at_y(y.getPoint());
        }
        throw std::runtime_error("neither point has fraction 0.0");
    }

};