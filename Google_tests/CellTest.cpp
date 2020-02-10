#include "gtest/gtest.h"
#include "../src/IntegralFrechet/Cell.h"
#include "../src/IntegralFrechet/metrics/include.h"

void AssertPointsEqual(const Point& a, const Point& b) {
    ASSERT_PRED3(approx_equal<Point>, a, b, ABS_TOL);
}

void AssertPathsEqual(const Points& a, const Points& b) {
    ASSERT_EQ(a.size(), b.size());
    for (size_t i = 0; i < a.size(); ++i) {
        AssertPointsEqual(a[i], b[i]);
    }
}

TEST(ComputeMatchingTest, LInfinity) {
    Cell cell({}, {}, {}, {});
    Points matching;

    cell = Cell({0, 0}, {0, 1}, {1, 1}, {1, 0});
    matching = compute_matching<ParamMetric::LInfinity_NoShortcuts>(cell, cell.s, cell.t);
    AssertPathsEqual(matching, {
        {0,                 0},
        {0.707106781186548, 0.707106781186548},
        {1.4142135623731,   1.4142135623731}
    });

    cell = Cell({0, 0}, {0, 0.5}, {1, 1}, {1, 0.5});
    matching = compute_matching<ParamMetric::LInfinity_NoShortcuts>(cell, cell.s, cell.t);
    AssertPathsEqual(matching,{
        {0, 0},
        {0.707106781186548, 0.5},
        {1.4142135623731, 1}
    });

    cell = Cell({0, 0}, {0, 0}, {1, 1}, {0, 3});
    matching = compute_matching<ParamMetric::LInfinity_NoShortcuts>(cell, cell.s, cell.t);
    AssertPathsEqual(matching, {
        {0, 0},
        {1.4142135623731, 2},
        {1.4142135623731, 3},
    });
}