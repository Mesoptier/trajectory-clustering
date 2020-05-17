#include "gtest/gtest.h"
#include "../src/IntegralFrechet/Cell.h"

TEST(CellTest, EllipseShape) {
    {
        // Parallel, same direction
        const Cell cell(
            {0, 0}, {0, 0},
            {1, 0}, {1, 0}
        );
        ASSERT_DOUBLE_EQ(cell.l, -1);
    }
    {
        // Parallel, opposite direction
        const Cell cell(
            {0, 0}, {0, 0},
            {1, 0}, {-1, 0}
        );
        ASSERT_DOUBLE_EQ(cell.l, 1);
    }
    {
        // Orthogonal
        const Cell cell(
            {0, 0}, {0, 0},
            {1, 0}, {0, 1}
        );
        ASSERT_DOUBLE_EQ(cell.l, 0);
    }
    {
        // Somewhere in between
        const Cell cell(
            {0, 0}, {0, 0},
            {1, 0}, {1, 1}
        );
        ASSERT_DOUBLE_EQ(cell.l, std::sqrt(2) / -2);
    }
}