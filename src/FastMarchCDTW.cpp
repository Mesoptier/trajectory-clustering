#include "FastMarchCDTW.h"
#include <armadillo>
#include <boost/heap/fibonacci_heap.hpp>

using namespace arma;

double
FastMarchCDTW::compute(
    const Curve<double>& curve1,
    const Curve<double>& curve2,
    double h,
    int imageNorm,
    int paramNorm,
    bool saveMatrices
) {
    double hi = h;
    double hj = h;
    unsigned int n_rows = curve1.getLength() / hi;
    unsigned int n_cols = curve2.getLength() / hj;

    // Mesh of the f-function
    dmat f_mat(n_rows, n_cols);

    // TODO: calculate this in a more clever way, should be possible to do linear in grid size
    for (unsigned int i = 0; i < n_rows; ++i) {
        for (unsigned int j = 0; j < n_cols; ++j) {
            f_mat(i, j) = norm(
                curve1.interp((double) i / (n_rows - 1)) - curve2.interp((double) j / (n_cols - 1)),
                imageNorm
            );
        }
    }

    Mat<short> tags(n_rows, n_cols);
    tags.fill(Far);

    dmat u_mat(n_rows, n_cols);
    u_mat.fill(INFINITY);

    // Fibonacci heap for efficient retrieval of considered point with lowest cost
    typedef boost::heap::fibonacci_heap<Node, boost::heap::compare<CompareNode>> Heap;
    Heap considered_points;
    std::map<Point, Heap::handle_type> handles;

    const Point start = {0, 0};
    const double start_cost = 0;
    tags(start.first, start.second) = Considered;
    u_mat(start.first, start.second) = start_cost;
    handles.emplace(start, considered_points.emplace(start, start_cost));

    const std::vector<Point> offsets = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};

    do {
        auto trial_node = considered_points.top();
        auto trial = trial_node.point;
        considered_points.pop();
        handles.erase(trial_node.point);
        tags(trial.first, trial.second) = Accepted;

        // Update direct neighbors
        for (const auto& offset : offsets) {
            const auto neighbor = trial + offset;
            if (!inBounds(neighbor, n_rows, n_cols)) {
                continue;
            }

            int i = neighbor.first;
            int j = neighbor.second;

            auto& tag = tags(i, j);
            if (tag == Accepted) {
                continue;
            }

            // Compute updated U value
            double u = eikonalUpdate(u_mat, f_mat, i, j, hi, hj, paramNorm);

            const double prev_u = u_mat(i, j);
            if (u < prev_u) {
                u_mat(i, j) = u;
            } else {
                u = prev_u;
            }

            if (tag == Far) {
                tag = Considered;
                handles.emplace(neighbor, considered_points.emplace(neighbor, u));
            } else if (tag == Considered && u < prev_u) {
                auto handle = handles[neighbor];
                considered_points.update(handle, Node{neighbor, u});
            }
        }
    } while (!considered_points.empty());

    if (saveMatrices) {
        // Compute gradient in I and J directions
        mat grad_i(n_rows, n_cols);
        mat grad_j(n_rows, n_cols);

        for (int i = 0; i < n_rows; ++i) {
            for (int j = 0; j < n_cols; ++j) {
                if (i > 0) {
                    grad_i(i, j) = (u_mat(i, j) - u_mat(i - 1, j)) / hi;
                }
                if (j > 0) {
                    grad_j(i, j) = (u_mat(i, j) - u_mat(i, j - 1)) / hj;
                }
            }
        }

        // Mix of gradient/coordinate descent to recover shortest path
        double step_size = 0.05;
        double precision = step_size;
        int max_iters = 100;

        rowvec current_x;
        rowvec next_x = {curve1.getLength(), curve2.getLength()};
        rowvec target_x = {0, 0};

        vec domain_i = linspace(0, curve1.getLength(), n_rows);
        vec domain_j = linspace(0, curve2.getLength(), n_cols);

        mat dx_i;
        mat dx_j;

        std::vector<rowvec> path_list;

        for (int i = 0; i < max_iters; ++i) {
            current_x = next_x;
            path_list.push_back(current_x);

            if (norm(target_x - current_x, paramNorm) < precision) {
                break;
            }

            interp2(domain_j, domain_i, grad_i, current_x.col(1), current_x.col(0), dx_i);
            interp2(domain_j, domain_i, grad_j, current_x.col(1), current_x.col(0), dx_j);

            const std::vector<rowvec> trial_offsets = {
                {-1, 0},
                {-dx_i[0], -dx_j[0]},
                {0, -1}
            };

            double min_u = INFINITY;
            for (const auto& offset : trial_offsets) {
                rowvec trial_x = current_x + step_size * normalise(offset, paramNorm);
                mat trial_u;
                interp2(domain_j, domain_i, u_mat, trial_x.col(1), trial_x.col(0), trial_u, "linear", INFINITY);

                if (trial_u[0] < min_u) {
                    min_u = trial_u[0];
                    next_x = trial_x;
                }
            }
        }

        // Create matching matrix from list of path vertices
        mat matching(path_list.size(), 2);
        for (int i = 0; i < path_list.size(); ++i) {
            matching.row(path_list.size() - i - 1) = path_list[i];
        }

        // Save the resulting matrices
        f_mat.save("f_mat.csv", csv_ascii);
        u_mat.save("u_mat.csv", csv_ascii);
        matching.save("matching.csv", csv_ascii);

        curve1.getVertices().save("curve1.csv", csv_ascii);
        curve2.getVertices().save("curve2.csv", csv_ascii);
    }

    return u_mat(n_rows - 1, n_cols - 1);
}

bool FastMarchCDTW::inBounds(FastMarchCDTW::Point point, unsigned int n_rows, unsigned int n_cols) {
    return (point.first >= 0) && (point.first < n_rows) && (point.second >= 0) && (point.second < n_cols);
}

double FastMarchCDTW::eikonalUpdate(const mat& u_mat, const mat& f_mat, int i, int j, double hi, double hj, int norm) {
    const double f = f_mat(i, j);

    if (f <= 0) {
        // TODO: Allow f = 0?
        throw std::range_error("f must be positive at all mesh vertices");
    }

    // Get u at neighbors
    double ui = INFINITY;
    if (i > 0) {
        ui = std::min(ui, u_mat(i - 1, j));
    }
    if (i < u_mat.n_rows - 1) {
        ui = std::min(ui, u_mat(i + 1, j));
    }

    double uj = INFINITY;
    if (j > 0) {
        uj = std::min(uj, u_mat(i, j - 1));
    }
    if (j < u_mat.n_cols - 1) {
        uj = std::min(uj, u_mat(i, j + 1));
    }

    // L1 norm
    if (norm == 1) {
        return std::min(ui + f * hi, uj + f * hj);
    }

    // L2 norm
    // TODO: This does not take into account distance travelled along curve
    if (norm == 2) {
        // Lower-dimensional update
        if (f < (uj - ui) / (hi * hi) && ui < uj) {
            return f * hi * hi + f * hi * hj + ui;
        } else if (f < (ui - uj) / (hj * hj) && uj < ui) {
            return f * hi * hj + f * hj * hj + uj;
        }

        // Solved in Mathematica: Solve[(Max[u-ui,0]/hi)^2+(Max[u-uj,0]/hj)^2==f^2&&hi>0&&hj>0&&f>0,u,Reals]
        const double hi2 = hi * hi;
        const double hi4 = hi2 * hi2;
        const double hj2 = hj * hj;
        const double hj4 = hj2 * hj2;
        const double f2 = f * f;
        const double ui2 = ui * ui;
        const double uj2 = uj * uj;

        return (hj2 * ui + hi2 * uj) / (hi2 + hj2) + sqrt(
            (f2 * hi4 * hj2 + f2 * hi2 * hj4 - hi2 * hj2 * ui2 + 2 * hi2 * hj2 * ui * uj - hi2 * hj2 * uj2)
            / ((hi2 + hj2) * (hi2 + hj2))
        );
    }

    throw std::logic_error("unsupported norm");
}
