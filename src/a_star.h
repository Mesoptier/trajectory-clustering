#pragma once

#include <queue>
#include <ostream>
#include <unordered_map>
#include <algorithm>
#include <iostream>
#include "io.h"

//#define A_STAR_LOGGING
#define A_STAR_STATS

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
        Points nodes_as_points;

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

    #ifdef A_STAR_STATS
    a_star::Stats stats{};
    #endif

    using QueueNode = std::pair<cost_t, Node>;
    std::priority_queue<QueueNode, std::vector<QueueNode>, std::greater<>> open_set;
    open_set.emplace(0, start);

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
            #endif
            #ifdef A_STAR_STATS
            std::cout << stats << " open_set remaining: " << open_set.size() << std::endl;
            io::export_points("data/out/debug_points.csv", stats.nodes_as_points);
            #endif
            return {g_score[goal], reconstruct_path<Graph>(came_from, goal)};
        }

        if (get_with_default(f_score, current, inf) < current_f) {
            // Already handled this node with a lower f score
            #ifdef A_STAR_LOGGING
            std::cout << " -> continue\n";
            #endif
            #ifdef A_STAR_STATS
            stats.nodes_skipped++;
            #endif
            continue;
        }

        #ifdef A_STAR_STATS
        stats.nodes_handled++;
        stats.nodes_as_points.push_back(graph.node_as_point(current));
        #endif
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

                #ifdef A_STAR_STATS
                stats.nodes_opened++;
                #endif
            }
        }
        neighbors.clear();
    }

    std::cout << "[A*] failure\n";
    throw std::runtime_error("A* failed to find a path");
}

template<class Graph>
std::pair<typename Graph::cost_t, std::vector<typename Graph::Node>>
bidirectional_dijkstra_search(const Graph& graph, typename Graph::Node s, typename Graph::Node t) {
    using Node = typename Graph::Node;
    using cost_t = typename Graph::cost_t;
    constexpr cost_t inf = std::numeric_limits<cost_t>::infinity();

    #ifdef A_STAR_STATS
    a_star::Stats stats{};
    #endif

    using QueueNode = std::pair<cost_t, Node>;
    std::priority_queue<QueueNode, std::vector<QueueNode>, std::greater<>> open_set_f;
    std::priority_queue<QueueNode, std::vector<QueueNode>, std::greater<>> open_set_b;
    open_set_f.emplace(0, s);
    open_set_b.emplace(0, t);

    std::unordered_map<Node, cost_t> cost_f;
    std::unordered_map<Node, cost_t> cost_b;
    cost_f[s] = 0;
    cost_b[t] = 0;

    std::unordered_map<Node, Node> came_from_f;
    std::unordered_map<Node, Node> came_from_b;

    cost_t lowest_cost = inf;
    Node lowest_cost_node = s;

    std::vector<Node> neighbors;

    BFDirection dir = BFDirection::Forward;

    while (!open_set_f.empty() && !open_set_b.empty()) {
        auto top_f = open_set_f.top();
        auto top_b = open_set_b.top();

        if (top_f.first + top_b.first >= lowest_cost) {
            #ifdef A_STAR_STATS
            std::cout << stats << " open_set remaining: " << (open_set_f.size() + open_set_b.size()) << std::endl;
            io::export_points("data/out/debug_points.csv", stats.nodes_as_points);
            #endif

            auto path_f = reconstruct_path<Graph>(came_from_f, lowest_cost_node);
            auto path_b = reconstruct_path<Graph>(came_from_b, lowest_cost_node);
            path_f.insert(path_f.end(), path_b.rbegin() + 1, path_b.rend());
            return {lowest_cost, path_f};
        }

//        dir = (open_set_f.size() > open_set_b.size()) ? BFDirection::Backward : BFDirection::Forward;
        dir = (dir == BFDirection::Forward) ? BFDirection::Backward : BFDirection::Forward;

        auto& open_set = dir == BFDirection::Forward ? open_set_f : open_set_b;
        auto& cost = dir == BFDirection::Forward ? cost_f : cost_b;
        auto& opposite_cost = dir == BFDirection::Forward ? cost_b : cost_f;
        auto& came_from = dir == BFDirection::Forward ? came_from_f : came_from_b;

        auto [current_cost, current] = open_set.top();
        open_set.pop();

        if (get_with_default(cost, current, inf) < current_cost) {
            #ifdef A_STAR_STATS
            stats.nodes_skipped++;
            #endif
            continue;
        }

        #ifdef A_STAR_STATS
        stats.nodes_handled++;
        stats.nodes_as_points.push_back(graph.node_as_point(current));
        #endif

        graph.get_neighbors(current, neighbors, dir);
        for (const Node& neighbor : neighbors) {
            cost_t neighbor_cost = dir == BFDirection::Forward
                ? current_cost + graph.cost(current, neighbor)
                : current_cost + graph.cost(neighbor, current);
            if (neighbor_cost < get_with_default(cost, neighbor, inf)) {
                came_from[neighbor] = current;
                cost[neighbor] = neighbor_cost;
                open_set.emplace(neighbor_cost, neighbor);

                #ifdef A_STAR_STATS
                stats.nodes_opened++;
                #endif

                cost_t path_cost = neighbor_cost + get_with_default(opposite_cost, neighbor, inf);
                if (path_cost < lowest_cost) {
                    lowest_cost = path_cost;
                    lowest_cost_node = neighbor;
                }
            }
        }
        neighbors.clear();
    }

    throw std::runtime_error("failed to find a path");
}