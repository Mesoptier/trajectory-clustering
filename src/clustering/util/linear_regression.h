#pragma once

#include "../../geom.h"
#include <Eigen/Dense>

using Eigen::MatrixXd;

Line linear_regression(Points& points, std::vector<distance_t>& weights) {
    size_t n = points.size();

    MatrixXd X(n, 2);
    MatrixXd y(n, 1);
    MatrixXd W(n, n);

    for (int i = 0; i < n; ++i) {
        X(i, 0) = 1;
        X(i, 1) = points[i].x;
        y(i, 0) = points[i].y;
        for (int j = 0; j < n; ++j) {
            if (i == j)
                W(i, j) = weights[i];
            else
                W(i, j) = 0;
        }
    }

    MatrixXd W_inv = W.inverse();
    MatrixXd X_trans = X.transpose();

    MatrixXd beta = (X_trans * W_inv * X).inverse() * X_trans * W_inv * y;

    if (std::isnan(beta(1, 0))) {
        std::cout << W << "\n";
        std::cout << "this is a problem...\n";
    }

    return Line::fromPointAndSlope({0, beta(0, 0)}, beta(1, 0));
}