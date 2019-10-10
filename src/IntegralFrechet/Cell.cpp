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
}

template
class Cell<double>;