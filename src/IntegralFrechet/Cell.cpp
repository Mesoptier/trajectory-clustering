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

    // Analyze edges
    bool parallel = approx_equal(edge1.slope, edge2.slope);

    // Monotone ellipse axis
    V slope = 1;
    V in1_intercept;

    if (!parallel) {
        // 1. find intersection between two lines
        // 2. convert intersection point to point in parameter space
        // 3. find intercept with in1/in2/out1/out2 axes
    } else {
        // 1. find closest point on edge2 to edge1.first
        // 2. convert closest point to point in parameter space
        // 3. -> automatically found intercept with in2 axis (because edge1.first equals 0 along in1 axis)

//        const auto norm_diff = arma::normalise(edge2.diff, 2);
//        in1_intercept = edge2.first + norm_diff * arma::dot(edge1.first - edge2.first, norm_diff);
    }

    // Default output values
    out1->fill(INFINITY);
    out2->fill(INFINITY);

    // Compute output values
    for (int i = 0; i < n1 + n2; ++i) {
        for (int o = 0; o < n1 + n2; ++o) {
            outValue(o) = std::min(outValue(o), inValue(i) + computeCost(i, o));
            // TODO: Keep track of in-point with shortest path to this out-point
        }
    }
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
V Cell<V>::computeCost(int i, int o) const {
    // 1. find shortest path (should be constant time, based on monotone ellipse axis)
    // 2. integrate over the 0-3 linear parts in the shortest path

    return integrate({0, 0}, {edge1.getLength(), edge1.getLength()});
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