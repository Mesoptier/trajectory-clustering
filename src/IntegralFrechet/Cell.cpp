#include "Cell.h"

template<class V>
Cell<V>::Cell(
    const Edge<V>& edge1,
    const Edge<V>& edge2,
    int n1,
    int n2,
    const std::shared_ptr<const arma::Col<V>>& in1,
    const std::shared_ptr<const arma::Col<V>>& in2
): edge1(edge1), edge2(edge2),
   n1(n1), n2(n2),
   in1(in1), in2(in2),
   out1(std::make_shared<arma::Col<V>>(n1)), out2(std::make_shared<arma::Col<V>>(n2)) {

    // === ANALYZE EDGES ===
    bool parallel = approx_equal(edge1.slope, edge2.slope);
    slope = 1;

    if (!parallel) {
        // 1. find intersection between two lines
        const auto imagePoint = intersectInfinite(edge1, edge2);
//        std::cout << imagePoint << std::endl;

        // 2. convert intersection point to point in parameter space
        midPoint = {edge1.param(imagePoint), edge2.param(imagePoint)};
//        std::cout << midPoint << std::endl;
    } else {
        // 1. find point on edge2 closest to edge1.first
        const auto norm_diff = arma::normalise(edge2.diff, 2);
        const auto imagePoint = edge2.first + norm_diff * arma::dot(edge1.first - edge2.first, norm_diff);

        // 2. convert closest point to point in parameter space
        midPoint = {0, edge2.param(imagePoint)};
    }

    // 3. find intercept
    intercept = -slope * midPoint(0) + midPoint(1);

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
            outValue(o) = std::min(outValue(o), inValue(i) + cost);
            std::cout << "inValue=" << inValue(i) << " cost=" << cost << " outValue=" << outValue(o) << std::endl;
            // TODO: Keep track of in-point with shortest path to this out-point
        }
    }
}

template<class V>
V Cell<V>::getResult() const {
    return (*out2)(n2 - 1);
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
V Cell<V>::computeCost(int i, int o) const {
    const auto a = inPoint(i);
    const auto b = outPoint(o);

    // TODO: Fix for loop so this case doesn't happen
    if (!(a(0) <= b(0) && a(1) <= b(1))) {
        return INFINITY;
    }

    // 1. find shortest path (should be constant time, based on monotone ellipse axis)
    Vertex<V> c1, c2;

    if (a(1) < a(0) * slope + intercept) {
        // Vertical
        c1 = {a(0), std::min(a(0) * slope + intercept, b(1))};
    } else if (a(1) > a(0) * slope + intercept) {
        // Horizontal
        c1 = {std::min((a(1) - intercept) / slope, b(0)), a(1)};
    } else {
        // On monotone ellipse axis
        c1 = a;
    }

    if (b(1) > b(0) * slope + intercept) {
        // Vertical
        c2 = {b(0), std::max(b(0) * slope + intercept, c1(1))};
    } else if (b(1) < b(0) * slope + intercept) {
        // Horizontal
        c2 = {std::max((b(1) - intercept) / slope, c1(0)), b(1)};
    } else {
        // On monotone ellipse axis
        c2 = b;
    }

    std::cout << std::endl << a << c1 << c2 << b;

    // 2. integrate over the 0-3 linear parts in the shortest path
    return integrate(a, c1) + integrate(c1, c2) + integrate(c2, b);
}

template<class V>
V Cell<V>::integrate(arma::Row<V> p1, arma::Row<V> p2) const {
    const Vertex<V> d1 = edge1.interpLength(p1(0)) - edge2.interpLength(p1(1));
    const Vertex<V> d2 = edge1.interpLength(p2(0)) - edge2.interpLength(p2(1));

    const V a = pow(d1(0) - d2(0), 2) + pow(d1(1) - d2(1), 2);
    const V b = 2 * (d1(0) * d2(0) - d1(0) * d1(0) + d1(1) * d2(1) - d1(1) * d1(1));
    const V c = d1(0) * d1(0) + d1(1) * d1(1);

    const V dist = arma::norm(p2 - p1, 1);

//    std::cout << "a=" << a << " b=" << b << " c=" << c << " dist=" << dist << std::endl;

    if (a != 0) {
        if (b != 0) {
            const V sa = sqrt(a);
            const V sc = sqrt(c);
            const V sabc = sqrt(a + b + c);

            return (
                2 * sa * (2 * a * sabc + b * (-sc + sabc))
                + (b * b - 4 * a * c) * (log(b + 2 * sa * sc) - log(2 * a + b + 2 * sa * sabc))
            ) / (8 * pow(a, 1.5)) * dist;
        }

        if (c != 0) {
            const V sa = sqrt(a);
            const V sac = sqrt(a + c);

            return (
                (sac / 2)
                - (c * (log(a) + log(c) - 2 * log(a + sa * sac))) / (4 * sa)
            ) * dist;
        }

        return sqrt(a) / 2 * dist;
    }

    if (b != 0) {
        throw std::runtime_error("unhandled integration case");
//        return 0 * dist;
    }

    return sqrt(c) * dist;
}

template
class Cell<double>;