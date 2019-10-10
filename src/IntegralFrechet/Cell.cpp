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
double Cell<V>::computeCost(int i, int o) const {
    // TODO: Compute cost
    return 0;
}

template
class Cell<double>;