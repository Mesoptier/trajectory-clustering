#pragma once

#include "../../geom.h"
// #include <Eigen/Dense>

// using Eigen::MatrixXd;

Line linear_regression(Points& points, std::vector<distance_t>& weights) {
    // size_t n = points.size();

    // MatrixXd X(n, 2);
    // MatrixXd y(n, 1);
    // MatrixXd W(n, n);

    // for (int i = 0; i < n; ++i) {
    //     X(i, 0) = 1;
    //     X(i, 1) = points[i].x;
    //     y(i, 0) = points[i].y;
    //     for (int j = 0; j < n; ++j) {
    //         if (i == j)
    //             W(i, j) = weights[i];
    //         else
    //             W(i, j) = 0;
    //     }
    // }

    // MatrixXd W_inv = W.inverse();
    // MatrixXd X_trans = X.transpose();

    // MatrixXd beta = (X_trans * W_inv * X).inverse() * X_trans * W_inv * y;

    // return Line::fromPointAndSlope({0, beta(0, 0)}, beta(1, 0));
    return Line::fromTwoPoints({0, 0}, {1, 0});
}

Line linear_regression_3d(Points& points, std::vector<distance_t>& weights, std::vector<distance_t>& heights) {
    // size_t n = points.size();

    // MatrixXd X(n, 3);
    // MatrixXd z(n, 1);
    // MatrixXd W(n, n);

    // for (int i = 0; i < n; ++i) {
    //     X(i, 0) = 1;
    //     X(i, 1) = points[i].x;
    //     X(i, 2) = points[i].y;
    //     z(i, 0) = heights[i];
    //     for (int j = 0; j < n; ++j) {
    //         if (i == j)
    //             W(i, j) = weights[i];
    //         else
    //             W(i, j) = 0;
    //     }
    // }

    // MatrixXd W_inv = W.inverse();
    // MatrixXd X_trans = X.transpose();

    // MatrixXd beta = (X_trans * W_inv * X).inverse() * X_trans * W_inv * z;
    // return Line::fromPointAndSlope({0, -beta(0, 0) / beta(2, 0)}, -beta(1, 0) / beta(2, 0));
    return Line::fromTwoPoints({0, 0}, {1, 0});
}