#include "FastMarchIntegralFrechet.h"
#include <armadillo>
#include <boost/heap/fibonacci_heap.hpp>

using namespace arma;

FastMarchIntegralFrechet::FastMarchIntegralFrechet(
    const Curve<double>& curve1, const Curve<double>& curve2, double h, int imageNorm, int paramNorm
) :
    curve1(curve1),
    curve2(curve2),
    imageNorm(imageNorm),
    paramNorm(paramNorm),

    // Compute size of matrices
    n_rows(ceil(curve1.getLength() / h)),
    n_cols(ceil(curve2.getLength() / h)),

    // Compute adjusted uniform step sizes
    hi(curve1.getLength() / n_rows),
    hj(curve2.getLength() / n_cols),

    f_mat(n_rows, n_cols),
    u_mat(n_rows, n_cols) {}

double FastMarchIntegralFrechet::computeDistance() {
    // Fill f-matrix
    // TODO: calculate this in a more clever way, should be possible to do linear in grid size
    for (unsigned int i = 0; i < n_rows; ++i) {
        for (unsigned int j = 0; j < n_cols; ++j) {
            f_mat(i, j) = norm(
                curve1.interp((double) i / (n_rows - 1)) - curve2.interp((double) j / (n_cols - 1)),
                imageNorm
            );
        }
    }

    // Fill u-matrix
    u_mat.fill(INFINITY);

    // Initialize and fill tags-matrix
    Mat<short> tags(n_rows, n_cols);
    tags.fill(Far);

    // Fibonacci heap for efficient retrieval of considered point with lowest cost
    typedef boost::heap::fibonacci_heap<Node, boost::heap::compare<CompareNode>> Heap;
    Heap considered_points;
    std::map<Point, Heap::handle_type> handles;

    const Point start = {0, 0};
    const double start_cost = 0;
    tags(start.first, start.second) = Considered;
    u_mat(start.first, start.second) = start_cost;
    handles.emplace(start, considered_points.emplace(start, start_cost));

    const std::vector<Point> offsets = {
        {-1, 0},
        {1,  0},
        {0,  -1},
        {0,  1}
    };

    do {
        auto trial_node = considered_points.top();
        auto trial = trial_node.point;
        considered_points.pop();
        handles.erase(trial_node.point);
        tags(trial.first, trial.second) = Accepted;

        // Update direct neighbors
        for (const auto& offset : offsets) {
            const auto neighbor = trial + offset;
            if (!inBounds(neighbor)) {
                continue;
            }

            int i = neighbor.first;
            int j = neighbor.second;

            auto& tag = tags(i, j);
            if (tag == Accepted) {
                continue;
            }

            // Compute updated U value
            double u = eikonalUpdate(i, j);

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

    return u_mat(n_rows - 1, n_cols - 1);
}

void FastMarchIntegralFrechet::computeMatching(double stepSize, int maxIterations) {
    // TODO: Allow configuring precision
    double precision = stepSize;

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
    rowvec current_x;
    rowvec next_x = {curve1.getLength(), curve2.getLength()};
    rowvec target_x = {0, 0};

    vec domain_i = linspace(0, curve1.getLength(), n_rows);
    vec domain_j = linspace(0, curve2.getLength(), n_cols);

    mat dx_i;
    mat dx_j;

    std::vector<rowvec> path_list;

    for (int i = 0; i < maxIterations; ++i) {
        current_x = next_x;
        path_list.push_back(current_x);

        if (norm(target_x - current_x, paramNorm) < precision) {
            break;
        }

        interp2(domain_j, domain_i, grad_i, current_x.col(1), current_x.col(0), dx_i);
        interp2(domain_j, domain_i, grad_j, current_x.col(1), current_x.col(0), dx_j);

        // Directions to check for potential next steps
        std::vector<rowvec> trial_offsets = {
            {-1,       0},
            {-dx_i[0], -dx_j[0]},
            {0,        -1},
        };

        // Add more directions with controllable precision
//            for (double j = 0; j <= 1; j += 0.01) {
//                trial_offsets.push_back({-j, -(1 - j)});
//            }

        double min_u = INFINITY;
        for (const auto& offset : trial_offsets) {
            rowvec trial_x = current_x + stepSize * normalise(offset, paramNorm);
            mat trial_u;
            interp2(domain_j, domain_i, u_mat, trial_x.col(1), trial_x.col(0), trial_u, "linear", INFINITY);

            if (trial_u[0] < min_u) {
                min_u = trial_u[0];
                next_x = trial_x;
            }
        }
    }

    // Create matching matrix from list of path vertices
    matching = mat(path_list.size(), 2);
    for (int i = 0; i < path_list.size(); ++i) {
        matching.row(path_list.size() - i - 1) = path_list[i];
    }
}

void FastMarchIntegralFrechet::computeCenter(double ratio) {
    center = mat(matching.n_rows, 2);

    for (int i = 0; i < matching.n_rows; ++i) {
        rowvec p1 = curve1.interpLength(matching(i, 0));
        rowvec p2 = curve2.interpLength(matching(i, 1));
        center.row(i) = p1 * ratio + p2 * (1 - ratio);
    }
}

void FastMarchIntegralFrechet::save() {
    f_mat.save("f_mat.csv", csv_ascii);
    u_mat.save("u_mat.csv", csv_ascii);
    matching.save("matching.csv", csv_ascii);
    center.save("center.csv", csv_ascii);

    curve1.getVertices().save("curve1.csv", csv_ascii);
    curve2.getVertices().save("curve2.csv", csv_ascii);
}

bool FastMarchIntegralFrechet::inBounds(Point point) {
    return (point.first >= 0) && (point.first < n_rows) && (point.second >= 0) && (point.second < n_cols);
}

double FastMarchIntegralFrechet::eikonalUpdate(int i, int j) {
    const double f = f_mat(i, j);

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
    if (paramNorm == 1) {
        return std::min(ui + f * hi, uj + f * hj);
    }

    if (f <= 0) {
        // TODO: Allow f = 0?
        throw std::range_error("f must be positive at all mesh vertices");
    }

    // L2 norm
    // TODO: This does not take into account distance travelled along curve
    if (paramNorm == 2) {
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

const mat& FastMarchIntegralFrechet::getCenter() const {
    return center;
}
