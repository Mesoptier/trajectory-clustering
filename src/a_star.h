#pragma once

#include <queue>

namespace {
    // From: https://stackoverflow.com/a/26958878/1639600
    template<class Map>
    const typename Map::mapped_type& get_with_default(
        const Map& m,
        const typename Map::key_type& key,
        const typename Map::mapped_type& defval
    ) {
        typename Map::const_iterator it = m.find(key);
        if (it == m.end())
            return defval;
        return it->second;
    }
}

template<class Graph>
void a_star_search(Graph graph, typename Graph::Node start, typename Graph::Node goal) {
    using Node = typename Graph::Node;
    using cost_t = typename Graph::cost_t;
    constexpr cost_t inf = std::numeric_limits<cost_t>::infinity();

    using QueueNode = std::pair<cost_t, Node>;
    std::priority_queue<QueueNode, std::vector<QueueNode>, std::greater<>> open_set;
    open_set.emplace(0, start);

    std::map<Node, Node> came_from;

    // For node n, g_score[n] is the cost of cheapest path from start to n
    std::map<Node, cost_t> g_score;
    g_score[start] = 0;

    // For node n, f_score[n] := g_score[n] + graph.heuristic_cost(n)
    std::map<Node, cost_t> f_score;
    f_score[start] = 0;

    std::vector<Node> neighbors;

    while (!open_set.empty()) {
        Node current;
        cost_t current_f;
        std::tie(current_f, current) = open_set.top();
        open_set.pop();

        std::cout << "[A*] node=" << current << " f=" << current_f;

        if (current == goal) {
            std::cout << " -> goal\n";
            return;
        }

        if (get_with_default(f_score, current, inf) < current_f) {
            // Already handled this node with a lower f score
            std::cout << " -> continue\n";
            continue;
        }

        std::cout << " -> handle\n";

        graph.get_neighbors(current, neighbors);
        for (const Node& neighbor : neighbors) {
            cost_t neighbor_g = g_score[current] + graph.cost(current, neighbor);
            if (neighbor_g < get_with_default(g_score, neighbor, inf)) {
                // Found better path to node neighbor
                came_from[neighbor] = current;
                g_score[neighbor] = neighbor_g;
                cost_t neighbor_f = neighbor_g + graph.heuristic_cost(neighbor, goal);
                f_score[neighbor] = neighbor_f;
                open_set.emplace(neighbor_f, neighbor);
            }
        }
    }

    std::cout << "[A*] failure\n";
}