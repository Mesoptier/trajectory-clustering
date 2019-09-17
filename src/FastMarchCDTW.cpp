#include "FastMarchCDTW.h"
#include <armadillo>
#include <boost/heap/fibonacci_heap.hpp>
#include <unordered_map>

using namespace arma;

void FastMarchCDTW::compute() {
    int n_rows = 5;
    int n_cols = 5;

    Mat<short> tags(n_rows, n_cols);
    tags.fill(Far);

    dmat costs(n_rows, n_cols);
    costs.fill(std::numeric_limits<double>::infinity());

    // Fibonacci heap for efficient retrieval of considered point with lowest cost
    typedef boost::heap::fibonacci_heap<Node, boost::heap::compare<CompareNode>> Heap;
    Heap considered_points;
    std::map<Point, Heap::handle_type> handles;

    Point start(0, 0);
    double start_cost = 0;
    tags(start.first, start.second) = Considered;
    costs(start.first, start.second) = start_cost;
    handles.emplace(start, considered_points.emplace(start, start_cost));

    std::cout << tags << std::endl;
    std::cout << costs << std::endl;

    do {
        auto trial_node = considered_points.top();
        auto trial = trial_node.point;
        considered_points.pop();
        handles.erase(trial_node.point);
        tags(trial.first, trial.second) = Accepted;

        for (int xx = -1; xx <= 1; ++xx) {
            for (int yy = -1; yy <= 1; ++yy) {
                // Only check horizontal/vertical neighbors
                if ((xx == 0) == (yy == 0)) {
                    continue;
                }

                // Check bounds
                if (trial.first + xx < 0 || trial.first + xx >= n_rows || trial.second + yy < 0 || trial.second + yy >= n_cols) {
                    continue;
                }

                Point neighbor(trial.first + xx, trial.second + yy);

                auto& tag = tags(neighbor.first, neighbor.second);

                if (tag == Accepted) {
                    continue;
                }

                // TODO: Recompute the cost (https://en.wikipedia.org/wiki/Eikonal_equation#Numerical_approximation)
                double cost = 1.0;
                costs(neighbor.first, neighbor.second) = cost;

                if (tag == Far) {
                    tag = Considered;
                    handles.emplace(neighbor, considered_points.emplace(neighbor, cost));
                } else { // tag == Considered
                    auto handle = handles[neighbor];
                    considered_points.update(handle, Node{neighbor, cost});
                }
            }
        }

        std::cout << tags << std::endl;
        std::cout << costs << std::endl;
    } while (!considered_points.empty());
}
