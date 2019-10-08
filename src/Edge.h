#ifndef CODE_EDGE_H
#define CODE_EDGE_H

#include <utility>
#include <ostream>
#include "Vertex.h"

template<class V>
class Edge : public std::pair<Vertex<V>, Vertex<V>> {
    using std::pair<Vertex<V>, Vertex<V>>::pair;

public:
    V getLength() const {
        return arma::norm(this->first - this->second, 2);
    }

    friend std::ostream& operator<<(std::ostream& os, const Edge& edge) {
        os << "V1: ";
        edge.first.print(os);
        os << "V2: ";
        edge.second.print(os);
        return os;
    }
};

#endif //CODE_EDGE_H
