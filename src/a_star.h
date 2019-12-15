#pragma once

#include <queue>
#include <ostream>
#include <unordered_map>

//#define A_STAR_LOGGING

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

namespace a_star {
    struct Stats {
        size_t nodes_opened;
        size_t nodes_handled;
        size_t nodes_skipped;

        friend std::ostream& operator<<(std::ostream& out, const Stats& stats) {
            out << "nodes_opened: " << stats.nodes_opened
                << " nodes_handled: " << stats.nodes_handled
                << " nodes_skipped: " << stats.nodes_skipped;
            return out;
        }
    };
}

template<class Graph>
std::vector<typename Graph::Node>
reconstruct_path(const std::unordered_map<typename Graph::Node, typename Graph::Node>& came_from, typename Graph::Node current) {
    std::vector<typename Graph::Node> path{current};

    auto it = came_from.find(current);
    while (it != came_from.end()) {
        path.push_back(it->second);
        it = came_from.find(it->second);
    }

    std::reverse(path.begin(), path.end());
    return path;
}

template<class Graph>
std::pair<typename Graph::cost_t, std::vector<typename Graph::Node>>
a_star_search(const Graph& graph, typename Graph::Node start, typename Graph::Node goal) {
    using Node = typename Graph::Node;
    using cost_t = typename Graph::cost_t;
    constexpr cost_t inf = std::numeric_limits<cost_t>::infinity();

    a_star::Stats stats{};

    using QueueNode = std::pair<cost_t, Node>;
    std::priority_queue<QueueNode, std::vector<QueueNode>, std::greater<>> open_set;
    open_set.emplace(0, start);
    stats.nodes_opened++;

    std::unordered_map<Node, Node> came_from;

    // For node n, g_score[n] is the cost of cheapest path from start to n
    std::unordered_map<Node, cost_t> g_score;
    g_score[start] = 0;

    // For node n, f_score[n] := g_score[n] + graph.heuristic_cost(n)
    std::unordered_map<Node, cost_t> f_score;
    f_score[start] = 0;

    std::vector<Node> neighbors;

    while (!open_set.empty()) {
        Node current;
        cost_t current_f;
        std::tie(current_f, current) = open_set.top();
        open_set.pop();

        #ifdef A_STAR_LOGGING
        std::cout << "[A*] node=" << current << " f=" << current_f;
        #endif

        if (current == goal) {
            #ifdef A_STAR_LOGGING
            std::cout << " -> goal\n";
            std::cout << stats << std::endl;
            #endif
            return {g_score[goal], reconstruct_path<Graph>(came_from, goal)};
        }

        if (get_with_default(f_score, current, inf) < current_f) {
            // Already handled this node with a lower f score
            stats.nodes_skipped++;
            #ifdef A_STAR_LOGGING
            std::cout << " -> continue\n";
            #endif
            continue;
        }

        stats.nodes_handled++;
        #ifdef A_STAR_LOGGING
        std::cout << " -> handle\n";
        #endif

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
                stats.nodes_opened++;
            }
        }
        neighbors.clear();
    }

    std::cout << "[A*] failure\n";
    throw std::runtime_error("A* failed to find a path");
}