#include <memory>
#include "Solver.h"

Solver::Solver(
    const Curve& curve1, const Curve& curve2,
    double h, ImageMetric imageMetric, ParamMetric paramMetric
) : curve1(curve1), curve2(curve2),
    n1(curve1.getNoVertices() - 1), n2(curve2.getNoVertices() - 1),
    imageMetric(imageMetric), paramMetric(paramMetric)
{

    // Create cell grid
    cells.reserve(n1 * n2);

    std::shared_ptr<arma::Col<distance_t>> in1;
    std::shared_ptr<arma::Col<distance_t>> in2;

    for (int i1 = 0; i1 < n1; ++i1) {
        const auto edge1 = curve1.getEdge(i1);
        unsigned int m1 = ceil(edge1.length / h) + 1;
        in1 = std::make_shared<arma::Col<distance_t>>(m1);
        in1->fill(INFINITY);
        if (i1 == 0) {
            in1->at(0) = 0;
        }

        for (int i2 = 0; i2 < n2; ++i2) {
            const auto edge2 = curve2.getEdge(i2);
            unsigned int m2 = ceil(edge2.length / h) + 1;

            if (i1 == 0) {
                in2 = std::make_shared<arma::Col<distance_t>>(m2);
                in2->fill(INFINITY);
            } else {
                in2 = cells[(i1 - 1) * n2 + i2].out2;
            }

            const Point offset = {curve1.getLength(i1), curve2.getLength(i2)};
            const Cell cell(edge1, edge2, m1, m2, in1, in2, offset, imageMetric, paramMetric);
            cells.push_back(cell);

            in1 = cell.out1;

            std::cout << cell << std::endl;
        }
    }
}

distance_t Solver::getDistance() const {
    return cells[n1 * n2 - 1].getResult();
}

Points Solver::getMatching() const {
    Point start = {curve1.getLength(), curve2.getLength()};
    Points matching = {start};

    int i1 = n1 - 1;
    int i2 = n2 - 1;
    bool done = false;

    while (!done) {
        // Do a final loop for the final cell
        done = i1 == 0 && i2 == 0;

        // Get cell
        auto cell = cells[i1 * n2 + i2];

        // Get cell offset
        const Point offset = cell.getOffset();

        // Get minimal path to the target point
        auto cellMatching = cell.getMinPath(start - offset);

        // Find index of next cell
        if (cellMatching[0].x == 0 && i1 > 0) {
            i1--;
        }
        if (cellMatching[0].y == 0 && i2 > 0) {
            i2--;
        }

        // Undo translation
        for (auto& point : cellMatching) {
            point += offset;
        }

        // Get current endpoint of the matching
        start = cellMatching.front();

        // Prepend minimal path to matching
        matching.insert(matching.end(), cellMatching.rbegin() + 1, cellMatching.rend());
    }

    std::reverse(matching.begin(), matching.end());
    return matching;
}

arma::Mat<distance_t> Solver::getBoundaryCosts() const {
    arma::Mat<distance_t> boundaryCosts(0, 3);
    for (const auto& cell : cells) {
        boundaryCosts = arma::join_cols(boundaryCosts, cell.getBoundaryCosts());
    }
    return boundaryCosts;
}
