#pragma once

#include <Curve.h>
#include <geom.h>
#include <ostream>
#include <utility>

class IntegralFrechet
{
private:
    const Curve curve1;
    const Curve curve2;

public:

    IntegralFrechet(Curve curve1, Curve curve2);

    void findPath();

    struct Cell {
        // Edges (edge1 = {p1, p1 + 1}, edge2 = {p2, p2 + 1})
        PointID p1;
        PointID p2;

        // TODO: Contains details about the cell:
        //  - Edges
        //  - EllH/EllV/etc.
    };

    //
    // Requirements for A* algorithm:
    //

    using cost_t = distance_t;
    using Node = CPosition;
    void get_neighbors(const Node& node, std::vector<Node>& neighbors) const;
    cost_t cost(const Node& s, const Node& t) const;
    cost_t heuristic_cost(const Node& s, const Node& goal) const;
};