#pragma once

#include <Curve.h>
#include <geom.h>
#include <ostream>
#include <utility>

class IntegralFrechet
{
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

    class Graph {
    private:
        const Curve curve1;
        const Curve curve2;

        Graph(Curve curve1, Curve curve2) : curve1(std::move(curve1)), curve2(std::move(curve2)) {}

    public:
        friend IntegralFrechet;

        using cost_t = distance_t;

        struct Node {
            // TODO: Replace with CPoint (currently only represents bottom-left corner of cell)
            PointID p1;
            PointID p2;

            Node() = default;
            Node(const PointID& p1, const PointID& p2) : p1(p1), p2(p2) {}

            bool operator==(const Node& rhs) const {
                return p1 == rhs.p1 && p2 == rhs.p2;
            }
            bool operator<(const Node& rhs) const {
                if (p1 < rhs.p1)
                    return true;
                if (rhs.p1 < p1)
                    return false;
                return p2 < rhs.p2;
            }
            bool operator>(const Node& rhs) const {
                return rhs < *this;
            }

            friend std::ostream& operator<<(std::ostream& os, const Node& node) {
                os << "(" << node.p1 << ", " << node.p2 << ")";
                return os;
            }
        };

        void get_neighbors(const Node& node, std::vector<Node>& neighbors) {
            if (node.p1 < curve1.size() - 1) {
                neighbors.emplace_back(node.p1 + 1, node.p2);

                if (node.p2 < curve2.size() - 1) {
                    neighbors.emplace_back(node.p1 + 1, node.p2 + 1);
                }
            }
            if (node.p2 < curve2.size() - 1) {
                neighbors.emplace_back(node.p1, node.p2 + 1);
            }
        }

        cost_t cost(const Node& s, const Node& t) {
            // TODO: Compute cost from s to t, assuming they are part of the same cell
            return 0;
        }

        cost_t heuristic_cost(const Node& s, const Node& t) {
            return 0;
        }
    };
};