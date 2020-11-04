#pragma once

#include <queue>
#include <ostream>
#include <unordered_map>
#include <algorithm>
#include <iostream>
#include "utils/io.h"

//#define A_STAR_LOGGING
#define A_STAR_STATS

namespace {
    // From: https://stackoverflow.com/a/26958878/1639600
    template<class Map>
    typename Map::mapped_type const& get_with_default(
        Map const& m,
        typename Map::key_type const& key,
        typename Map::mapped_type const& defval
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

        friend std::ostream& operator<<(std::ostream& out, Stats const& stats) {
            out << "nodes_opened: " << stats.nodes_opened
                << " nodes_handled: " << stats.nodes_handled
                << " nodes_skipped: " << stats.nodes_skipped;
            return out;
        }
    };
}

// namespace shortest_path_algs {

struct SearchStat {
    size_t nodes_opened;
    size_t nodes_handled;
    size_t nodes_skipped;
    Points nodes_as_points;
};

template<class Graph>
struct SearchResult {
    typename Graph::cost_t cost;
    std::vector<typename Graph::Node> path;
    SearchStat stat;
};

// }

// using namespace shortest_path_algs;

template<class Node>
std::vector<Node>
reconstruct_path(std::unordered_map<Node, Node> const& came_from, Node current) {
    std::vector<Node> path{current};

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
a_star_search(Graph const& graph, typename Graph::Node start, typename Graph::Node goal) {
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
        for (Node const& neighbor : neighbors) {
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
SearchResult<Graph>
bidirectional_dijkstra_search(Graph const& graph, typename Graph::Node s, typename Graph::Node t) {
    using Node = typename Graph::Node;
    using NodeID = size_t;
    using cost_t = typename Graph::cost_t;
    constexpr cost_t inf = std::numeric_limits<cost_t>::infinity();

    SearchStat stat{};

    std::vector<Node> nodes;
    std::unordered_map<Node, NodeID> node_ids;

    using QueueNode = std::pair<cost_t, NodeID>;
    std::priority_queue<QueueNode, std::vector<QueueNode>, std::greater<>> open_set_f;
    std::priority_queue<QueueNode, std::vector<QueueNode>, std::greater<>> open_set_b;

    std::unordered_map<NodeID, cost_t> cost_f;
    std::unordered_map<NodeID, cost_t> cost_b;

    std::unordered_map<NodeID, NodeID> came_from_f;
    std::unordered_map<NodeID, NodeID> came_from_b;

    NodeID s_id = nodes.size();
    nodes.push_back(s);
    node_ids.emplace(s, s_id);
    open_set_f.emplace(0, s_id);
    cost_f[s_id] = 0;

    NodeID t_id = nodes.size();
    nodes.push_back(t);
    node_ids.emplace(t, t_id);
    open_set_b.emplace(0, t_id);
    cost_b[t_id] = 0;

    cost_t lowest_cost = inf;
    NodeID lowest_cost_node_id = s_id;

    std::vector<Node> neighbors;
    std::vector<cost_t> costs;

    BFDirection dir = BFDirection::Forward;

    while (!open_set_f.empty() && !open_set_b.empty()) {
        auto top_f_cost = open_set_f.top().first;
        auto top_b_cost = open_set_b.top().first;

        if (top_f_cost + top_b_cost >= lowest_cost) {
            auto id_path_f = reconstruct_path<NodeID>(came_from_f, lowest_cost_node_id);
            auto id_path_b = reconstruct_path<NodeID>(came_from_b, lowest_cost_node_id);
            id_path_f.insert(id_path_f.end(), id_path_b.rbegin() + 1, id_path_b.rend());

            std::vector<Node> path;
            for (auto id : id_path_f) {
                path.push_back(nodes[id]);
            }

            return {
                lowest_cost,
                path,
                stat,
            };
        }

        dir = (open_set_f.size() > open_set_b.size()) ? BFDirection::Backward : BFDirection::Forward;
//        dir = (dir == BFDirection::Forward) ? BFDirection::Backward : BFDirection::Forward;

        auto& open_set = dir == BFDirection::Forward ? open_set_f : open_set_b;
        auto& cost = dir == BFDirection::Forward ? cost_f : cost_b;
        auto& opposite_cost = dir == BFDirection::Forward ? cost_b : cost_f;
        auto& came_from = dir == BFDirection::Forward ? came_from_f : came_from_b;

        auto [current_cost, current_id] = open_set.top();
        open_set.pop();
        auto current = nodes[current_id];

        if (get_with_default(cost, current_id, inf) < current_cost) {
            ++stat.nodes_skipped;
            continue;
        }
        ++stat.nodes_handled;

        #ifdef A_STAR_STATS
        stat.nodes_as_points.push_back(graph.node_as_point(current));
        #endif

        graph.get_neighbors(current, neighbors, costs, dir);
        for (size_t i = 0; i < neighbors.size(); ++i) {
            Node const& neighbor = neighbors[i];
            cost_t neighbor_cost = current_cost + costs[i];

            auto neighbor_id_it = node_ids.find(neighbor);
            NodeID neighbor_id = neighbor_id_it != node_ids.end() ? neighbor_id_it->second : nodes.size();

            if (neighbor_id == nodes.size()) {
                nodes.push_back(neighbor);
                node_ids[neighbor] = neighbor_id;
            }

            if (neighbor_cost < get_with_default(cost, neighbor_id, inf)) {
                came_from[neighbor_id] = current_id;
                cost[neighbor_id] = neighbor_cost;
                open_set.emplace(neighbor_cost, neighbor_id);

                ++stat.nodes_opened;

                cost_t path_cost = neighbor_cost + get_with_default(opposite_cost, neighbor_id, inf);
                if (path_cost < lowest_cost) {
                    lowest_cost = path_cost;
                    lowest_cost_node_id = neighbor_id;
                }
            }
        }
        neighbors.clear();
        costs.clear();
    }

    throw std::runtime_error("failed to find a path");
}
