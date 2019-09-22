#include "FastMarchCDTW.h"
#include <armadillo>
#include <boost/heap/fibonacci_heap.hpp>

using namespace arma;

double FastMarchCDTW::compute(const Curve<double>& curve1, const Curve<double>& curve2, double h, bool saveMatrices) {
    double hi = h;
    double hj = h;
    unsigned int n_rows = curve1.getLength() / hi;
    unsigned int n_cols = curve2.getLength() / hj;

    // Mesh of the f-function
    mat local_costs(n_rows, n_cols);

    // TODO: calculate this in a more clever way, should be possible to do linear in grid size
    for (unsigned int i = 0; i < n_rows; ++i) {
        for (unsigned int j = 0; j < n_cols; ++j) {
            local_costs(i, j) = norm(
                curve1.interp((double) i / (n_rows - 1)) - curve2.interp((double) j / (n_cols - 1)), 2);
        }
    }

    Mat<short> tags(n_rows, n_cols);
    tags.fill(Far);

    dmat costs(n_rows, n_cols);
    costs.fill(INFINITY);

    // Fibonacci heap for efficient retrieval of considered point with lowest cost
    typedef boost::heap::fibonacci_heap<Node, boost::heap::compare<CompareNode>> Heap;
    Heap considered_points;
    std::map<Point, Heap::handle_type> handles;

    const Point start = {0, 0};
    const double start_cost = 0;
    tags(start.first, start.second) = Considered;
    costs(start.first, start.second) = start_cost;
    handles.emplace(start, considered_points.emplace(start, start_cost));

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
                if (trial.first + xx < 0 || trial.first + xx >= n_rows || trial.second + yy < 0 ||
                    trial.second + yy >= n_cols) {
                    continue;
                }

                Point neighbor = {trial.first + xx, trial.second + yy};

                auto& tag = tags(neighbor.first, neighbor.second);

                if (tag == Accepted) {
                    continue;
                }

                // Solve in Mathematica: Solve[(Max[u-ui,0]/hi)^2+(Max[u-uj,0]/hj)^2==f^2&&hi>0&&hj>0&&f>0,u,Reals]
                double ui = INFINITY;
                if (neighbor.first > 0) {
                    ui = std::min(ui, costs(neighbor.first - 1, neighbor.second));
                }
                if (neighbor.first < n_rows - 1) {
                    ui = std::min(ui, costs(neighbor.first + 1, neighbor.second));
                }

                double uj = INFINITY;
                if (neighbor.second > 0) {
                    uj = std::min(uj, costs(neighbor.first, neighbor.second - 1));
                }
                if (neighbor.second < n_cols - 1) {
                    uj = std::min(uj, costs(neighbor.first, neighbor.second + 1));
                }

                const double f = local_costs(neighbor.first, neighbor.second);

                if (f <= 0) {
                    // TODO: Allow f = 0?
                    throw std::range_error("f must be positive at all mesh vertices");
                }

                double u;
                if (f < (uj - ui) / hi && ui < uj) {
                    u = f * hi + ui;
                } else if (f < (ui - uj) / hj && uj < ui) {
                    u = f * hj + uj;
                } else {
                    const double hi2 = hi * hi;
                    const double hi4 = hi2 * hi2;
                    const double hj2 = hj * hj;
                    const double hj4 = hj2 * hj2;
                    const double f2 = f * f;
                    const double ui2 = ui * ui;
                    const double uj2 = uj * uj;

                    u = (hj2 * ui + hi2 * uj) / (hi2 + hj2);
                    u += sqrt(
                        (f2 * hi4 * hj2 + f2 * hi2 * hj4 - hi2 * hj2 * ui2 + 2 * hi2 * hj2 * ui * uj - hi2 * hj2 * uj2)
                        / ((hi2 + hj2) * (hi2 + hj2))
                    );
                }

                const double prev_u = costs(neighbor.first, neighbor.second);
                if (u < prev_u) {
                    costs(neighbor.first, neighbor.second) = u;
                } else {
                    u = prev_u;
                }

                if (tag == Far) {
                    tag = Considered;
                    handles.emplace(neighbor, considered_points.emplace(neighbor, u));
                } else if (u < prev_u) { // tag == Considered
                    auto handle = handles[neighbor];
                    considered_points.update(handle, Node{neighbor, u});
                }
            }
        }
    } while (!considered_points.empty());

    // Compute gradient
    mat grad_i(n_rows, n_cols);
    mat grad_j(n_rows, n_cols);

    for (int i = 0; i < n_rows; ++i) {
        for (int j = 0; j < n_cols; ++j) {
            if (i > 0) {
                grad_i(i, j) = (costs(i, j) - costs(i - 1, j)) / hi;
            }
            if (j > 0) {
                grad_j(i, j) = (costs(i, j) - costs(i, j - 1)) / hj;
            }
        }
    }

    // Gradient descent to recover shortest path
    double step_size = std::min(hi, hj);
    double precision = 0.00001;
    int max_iters = 10000;

    rowvec current_x;
    rowvec next_x = {curve1.getLength(), curve2.getLength()};
    rowvec target_x = {0, 0};

    vec domain_i = linspace(0, curve1.getLength(), n_rows);
    vec domain_j = linspace(0, curve2.getLength(), n_cols);

    mat dx_i;
    mat dx_j;
    rowvec dx;

    std::vector<rowvec> path_list;

    for (int i = 0; i < max_iters; ++i) {
        current_x = next_x;
        path_list.push_back(current_x);

        if (norm(target_x - current_x) < precision) {
            break;
        }

        // Compute local gradient
        interp2(domain_j, domain_i, grad_i, current_x.col(1), current_x.col(0), dx_i);
        interp2(domain_j, domain_i, grad_j, current_x.col(1), current_x.col(0), dx_j);
        dx = {dx_i[0], dx_j[0]};

        next_x = current_x - step_size * normalise(dx, 2);

        next_x(0) = std::min(std::max(next_x(0), 0.0), curve1.getLength());
        next_x(1) = std::min(std::max(next_x(1), 0.0), curve2.getLength());
    }

    // Compute matching
    mat matching_param(path_list.size(), 2);
    mat matching_image(path_list.size(), 4);

    for (int i = 0; i < path_list.size(); ++i) {
        matching_param.row(i) = path_list[i];
        matching_image.row(i) = join_rows(curve1.interpLength(path_list[i](0)), curve2.interpLength(path_list[i](1)));
    }

    // Save the resulting matrices
    if (saveMatrices) {
        local_costs.save("mesh.csv", csv_ascii);
        costs.save("costs.csv", csv_ascii);
        matching_param.save("matching_param.csv", csv_ascii);
        matching_image.save("matching_image.csv", csv_ascii);

        curve1.getVertices().save("curve1.csv", csv_ascii);
        curve2.getVertices().save("curve2.csv", csv_ascii);
    }

    return costs(n_rows - 1, n_cols - 1);
}

bool FastMarchCDTW::inBounds(FastMarchCDTW::Point point, unsigned int n_rows, unsigned int n_cols) {
    return (point.first >= 0) && (point.first < n_rows) && (point.second >= 0) && (point.second < n_cols);
}
