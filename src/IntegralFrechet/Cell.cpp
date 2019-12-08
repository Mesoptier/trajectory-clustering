#include "Cell.h"

template<class V>
Cell<V>::Cell(
    const Edge<V>& edge1,
    const Edge<V>& edge2,
    int n1,
    int n2,
    const std::shared_ptr<const arma::Col<V>>& in1,
    const std::shared_ptr<const arma::Col<V>>& in2,
    Point offset,
    ImageMetric imageMetric,
    ParamMetric paramMetric
): edge1(edge1), edge2(edge2),
   n1(n1), n2(n2),
   in1(in1), in2(in2),
   offset(offset),
   imageMetric(imageMetric),
   paramMetric(paramMetric),
   out1(std::make_shared<arma::Col<V>>(n1)),
   out2(std::make_shared<arma::Col<V>>(n2)),
   outOrigin(n1 + n2) {

    // === ANALYZE EDGES ===
    if (!isParallel(edge1.line, edge2.line)) { // -> Infinite lines intersect somewhere
        // 1. find intersection between two lines
        const auto imagePoint = intersect(edge1.line, edge2.line);

        // 2. convert intersection point to point in parameter space
        midPoint = {edge1.param(imagePoint), edge2.param(imagePoint)};
    } else { // -> Infinite lines are parallel
        // TODO: Handle case where edges are in opposite directions

        // 1. find point on edge2 closest to edge1.first
        const Point imagePoint =
            edge2.first + edge2.line.direction * arma::dot(edge1.first - edge2.first, edge2.line.direction);

        // 2. convert closest point to point in parameter space
        midPoint = {0, edge2.param(imagePoint)};
    }

    // Initialise axes
    // TODO: Ellipse Axis might not have slope 1 if ImageMetric == L1...
    //       (Might be fixed by using L1 for all image-related distances?)
    ellipseAxis = Line::fromPointAndSlope(midPoint, 1);

    if (isPerpendicular(edge1.line, edge2.line)) {
        ellH = Line::fromPointAndSlope(midPoint, INFINITY);
        ellV = Line::fromPointAndSlope(midPoint, 0);
    } else {
        const auto dot = arma::dot(edge2.line.direction, edge1.line.direction);
        ellH = Line(midPoint, {dot, 1});
        ellV = Line(midPoint, {1, dot});
    }

    // === COMPUTE OUTPUT VALUES ===
    // Default output values
    out1->fill(INFINITY);
    out2->fill(INFINITY);

    // Compute output values
    for (int i = 0; i < n1 + n2; ++i) {
        for (int o = 0; o < n1 + n2; ++o) {
            if (std::isinf(inValue(i))) {
                continue;
            }

            const auto cost = computeCost(i, o);
            if (inValue(i) + cost < outValue(o)) {
                outValue(o) = inValue(i) + cost;
                outOrigin(o) = i;
            }
        }
    }
}

template<class V>
V Cell<V>::getResult() const {
    return (*out2)(n2 - 1);
}

template<class V>
arma::Mat<V> Cell<V>::getPath(int i, int o) const {
    const auto a = inPoint(i);
    const auto b = outPoint(o);

    if (paramMetric == ParamMetric::L1) {
        arma::Row<V> c1, c2;

        // NOTE: ellipseAxis is not vertical, so getY and getX are always defined

        if (ellipseAxis.includesPoint(a)) {
            // On monotone ellipse axis
            c1 = a;
        } else if (a(1) < ellipseAxis.getY(a(0))) {
            // Vertical
            c1 = {a(0), std::min(ellipseAxis.getY(a(0)), b(1))};
        } else {
            // Horizontal
            c1 = {std::min(ellipseAxis.getX(a(1)), b(0)), a(1)};
        }

        if (ellipseAxis.includesPoint(b)) {
            // On monotone ellipse axis
            c2 = b;
        } else if (b(1) > ellipseAxis.getY(b(0))) {
            // Vertical
            c2 = {b(0), std::max(ellipseAxis.getY(b(0)), c1(1))};
        } else {
            // Horizontal
            c2 = {std::max(ellipseAxis.getX(b(1)), c1(0)), b(1)};
        }

        arma::Mat<V> path(4, 2);
        path.row(0) = a;
        path.row(1) = c1;
        path.row(2) = c2;
        path.row(3) = b;
        return path;
    }

    if (paramMetric == ParamMetric::LInfinity_NoShortcuts) {
        // TODO: Follow steepest descent from a and b in monotone/reverse-monotone direction

    }

    throw std::logic_error("unsupported param metric");
}

template<class V>
Points Cell<V>::getMinPath(Point target) const {
    int targetIndex = outIndex(target);
    return getPath(outOrigin(targetIndex), targetIndex);
}

template<class V>
arma::Mat<V> Cell<V>::getBoundaryCosts() const {
    arma::Mat<V> boundaryCosts(n1 + n2, 3);
    for (int i = 0; i < n1 + n2; ++i) {
        const arma::Row<V> point = outPoint(i) + getOffset();
        boundaryCosts(i, 0) = point(0);
        boundaryCosts(i, 1) = point(1);
        boundaryCosts(i, 2) = outValue(i);
    }
    return boundaryCosts;
}

template<class V>
const V& Cell<V>::inValue(int i) const {
    return i < n1 ? (*in1)(i) : (*in2)(i - n1);
}

template<class V>
V& Cell<V>::outValue(int i) const {
    return i < n1 ? (*out1)(i) : (*out2)(i - n1);
}

template<class V>
arma::Row<V> Cell<V>::inPoint(int i) const {
    if (i < n1) {
        return {(V) i / (n1 - 1) * edge1.length, 0};
    } else {
        return {0, (V) (i - n1) / (n2 - 1) * edge2.length};
    }
}

template<class V>
arma::Row<V> Cell<V>::outPoint(int i) const {
    if (i < n1) {
        return {(V) i / (n1 - 1) * edge1.length, edge2.length};
    } else {
        return {edge1.length, (V) (i - n1) / (n2 - 1) * edge2.length};
    }
}

template<class V>
int Cell<V>::inIndex(arma::Row<V> p) const {
    return 0;
}

template<class V>
int Cell<V>::outIndex(arma::Row<V> p) const {
    if (approx_equal(p(1), edge2.length)) {
        // on out1 edge
        return round(p(0) * (n1 - 1) / edge1.length);
    } else {
        // on out2 edge
        return round(p(1) * (n2 - 1) / edge2.length) + n1;
    }
}

template<class V>
void Cell<V>::steepestDescent(PointsList& list, Point s, Point t) const {
    // TODO: Currently supports LInfinity_NoShortcuts only
    // TODO: Support degenerate cells (parallel edges)

    list.push_back(s);

    // Easy shortest path when the source and target are equal
    if (arma::approx_equal(s, t, "absdiff", ABS_TOL)) {
        return;
    }

    // If source is monotone greater than target, we need to do steepest
    // descent in reverse monotone direction.
    auto dir = MonotoneComparator::getDirection(s, t);
    MonotoneComparator compare(dir);

    // Multiply with line direction when checking on which side of the line a point lies
    int lineDir = (dir == MonotoneComparator::LowerFirst ? 1 : -1);

    // Boundaries of the subcell
    Line tHor(t, {1, 0});
    Line tVer(t, {0, 1});
    Line sHor(s, {1, 0});
    Line sVer(s, {0, 1});

    if (ellV.includesPoint(s) || ellH.includesPoint(s) || ellipseAxis.includesPoint(s)) {
        // -> On one of the steepest-descent diagonals
        if (compare(s, midPoint)) {
            // -> "Below" midPoint (= lowest point)

            // Which diagonal are we on again?
            Line line = ellV.includesPoint(s) ? ellV :
                        ellH.includesPoint(s) ? ellH :
                        ellipseAxis;

            list.push_back(std::min(
                {midPoint, intersect(line, tHor), intersect(line, tVer)},
                compare
            ));

            return;
        } else {
            // -> No steepest descent is possible
            return;
        }
    }

    // Above (left of) ellV
    if (perp(s - ellV.origin, lineDir * ellV.direction) < 0) {
        if (tVer.includesPoint(s)) {
            return;
        } else if (compare(s, midPoint)) {
            steepestDescent(list, std::min(intersect(ellV, sHor), intersect(tVer, sHor), compare), t);
            return;
        } else if (perp(s - ellH.origin, lineDir * ellH.direction) < 0) {
            steepestDescent(list, std::min(intersect(ellH, sHor), intersect(tVer, sHor), compare), t);
            return;
        } else {
            return;
        }
    }

    // Right of ellH
    if (perp(s - ellH.origin, lineDir * ellH.direction) > 0) {
        if (tHor.includesPoint(s)) {
            return;
        } else if (compare(s, midPoint)) {
            steepestDescent(list, std::min(intersect(ellH, sVer), intersect(tHor, sVer), compare), t);
            return;
        } else if (perp(s - ellV.origin, lineDir * ellV.direction) > 0) {
            steepestDescent(list, std::min(intersect(ellV, sVer), intersect(tHor, sVer), compare), t);
            return;
        } else {
            return;
        }
    }

    // Move parallel to ellipseAxis
    Line diag(s, ellipseAxis.direction);

    // Left of diagonal
    if (perp(s - ellipseAxis.origin, lineDir *ellipseAxis.direction) < 0) {
        steepestDescent(list, std::min(
            {intersect(ellV, diag), intersect(tHor, diag), intersect(tVer, diag)},
            compare
        ), t);
        return;
    } else {
        steepestDescent(list, std::min(
            {intersect(ellH, diag), intersect(tHor, diag), intersect(tVer, diag)},
            compare
        ), t);
        return;
    }
}

template<class V>
V Cell<V>::computeCost(int i, int o) const {
    const auto a = inPoint(i);
    const auto b = outPoint(o);

    // TODO: Fix for loop so this case doesn't happen
    if (!(a(0) <= b(0) && a(1) <= b(1))) {
        return INFINITY;
    }

    V cost = 0;
    const auto path = getPath(i, o);
    for (int j = 0; j < path.n_rows - 1; ++j) {
        cost += integrate(path.row(j), path.row(j + 1), imageMetric);
    }
    return cost;
}

template<class V>
V Cell<V>::integrate(arma::Row<V> p1, arma::Row<V> p2, ImageMetric imageMetric) const {
    // Get difference vectors between start- and endpoints of the two sub-edges.
    const Vertex<V> d1 = edge1.interpLength(p1(0)) - edge2.interpLength(p1(1));
    const Vertex<V> d2 = edge1.interpLength(p2(0)) - edge2.interpLength(p2(1));

    const V dx1 = d1(0);
    const V dy1 = d1(1);
    const V dx2 = d2(0);
    const V dy2 = d2(1);

    // Length of the edge in parameter space
    V dist;
    switch (paramMetric) {
        case ParamMetric::L1:
            dist = arma::norm(p2 - p1, 1);
            break;
        case ParamMetric::LInfinity_NoShortcuts:
            dist = arma::norm(p2 - p1, "inf");
            break;
    }

    return integrate(dx1, dy1, dx2, dy2, imageMetric) * dist;
}

template<class V>
V Cell<V>::integrate(V dx1, V dy1, V dx2, V dy2, ImageMetric imageMetric) const {
    if (imageMetric == ImageMetric::L1) {
        V cost = 0.0;

        // Mathematica: Integrate[Abs[(1-t)dx1+t dx2],{t,0,1}] (and similar for dy1 & dy2)
        if ((dx1 <= 0 && dx2 <= 0) || (dx1 >= 0 && dx2 >= 0)) {
            cost += (std::abs(dx1) + std::abs(dx2)) / 2;
        } else {
            cost += (pow(dx1, 2) + pow(dx2, 2)) / (2 * (std::abs(dx1) + std::abs(dx2)));
        }

        if ((dy1 <= 0 && dy2 <= 0) || (dy1 >= 0 && dy2 >= 0)) {
            cost += (std::abs(dy1) + std::abs(dy2)) / 2;
        } else {
            cost += (pow(dy1, 2) + pow(dy2, 2)) / (2 * (std::abs(dy1) + std::abs(dy2)));
        }

        return cost;
    }

    if (imageMetric == ImageMetric::L2) {
        if ((dx1 == 0 && dx2 == 0) || (dy1 == 0 && dy2 == 0)) {
            return integrate(dx1, dy1, dx2, dy2, ImageMetric::L1);
        }
    }

    // TODO: Clean up integration below this point

    // Coefficients are the same for both L2 and L2_Squared
    const V a = pow(dx1 - dx2, 2) + pow(dy1 - dy2, 2);
    const V b = 2 * (dx1 * dx2 - dx1 * dx1 + dy1 * dy2 - dy1 * dy1);
    const V c = dx1 * dx1 + dy1 * dy1;

    // === Image Metric: L2 Squared ===
    if (imageMetric == ImageMetric::L2_Squared) {
        return (a / 3 + b / 2 + c);
    }

    // else:
    // === Image Metric: L2 ===
    if (!approx_zero(a)) {
//        const V tmp1 = sqrt(a * (a + b + c));
//        const V tmp2 = sqrt(a * c);
//
//        return (
//            4 * a * tmp1 +
//                2 * b * (tmp1 - tmp2) +
//                (b * b - 4 * a * c) * (log(b + 2 * tmp2) - log(2 * a + b + 2 * tmp1))
//        ) / (8 * pow(a, 1.5));

        if (!approx_zero(b)) {
            if (b > 0) {
                const V sa = sqrt(a);
                const V sc = sqrt(c);
                const V sabc = sqrt(a + b + c);

                return (
                    2 * sa * (2 * a * sabc + b * (-sc + sabc))
                        + (b * b - 4 * a * c) * (log(b + 2 * sa * sc) - log(2 * a + b + 2 * sa * sabc))
                ) / (8 * pow(a, 1.5));
            } else {
                if (approx_equal(b, -c)) {
                    if (approx_equal(4 * a, c)) {
                        return (sqrt(a) + sqrt(c)) / 2;
                    }

                    return (
                        (4 * a * a)
                            - (2 * a * c)
                            + 2 * c * sqrt(a * c)
                            + c * (-4 * a + c) * (-log(4 * a - c) + log(b + 2 * sqrt(a * c)))
                    ) / (8 * pow(a, 1.5));
                }
            }
        }

        if (!approx_zero(c)) {
            const V sa = sqrt(a);
            const V sac = sqrt(a + c);

            return (
                (sac / 2)
                    - (c * (log(a) + log(c) - 2 * log(a + sa * sac))) / (4 * sa)
            );
        }

        return sqrt(a) / 2;
    }

    if (!approx_zero(b)) {
        return (
            2 * (-pow(c, 1.5) + pow(b + c, 1.5))
        ) / (3 * b);
    }

    return sqrt(c);
}

template
class Cell<double>;