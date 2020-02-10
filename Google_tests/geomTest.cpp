#include "gtest/gtest.h"
#include "../src/geom.h"

TEST(PointTest, Perp) {
    // Collinear:
    ASSERT_DOUBLE_EQ(perp({1, 1}, {5, 5}), 0.0);

    // Second vector to the left of first:
    ASSERT_GT(perp({1, 1}, {0, 2}), 0);

    // Second vector to the right of first:
    ASSERT_LT(perp({1, 1}, {2, 0}), 0);
}

TEST(MonotoneComparatorTest, GetDirection) {
    // Lower first:
    ASSERT_EQ(getMonotoneDirection({0, 0}, {1, 1}), BFDirection::Forward);
    // Higher first:
    ASSERT_EQ(getMonotoneDirection({1, 1}, {0, 0}), BFDirection::Backward);
    // Not monotone:
    ASSERT_ANY_THROW(getMonotoneDirection({1, 0}, {0, 1}));
    // Equal:
    ASSERT_ANY_THROW(getMonotoneDirection({0, 0}, {0, 0}));
}

TEST(MonotoneComparatorTest, Compare) {
    using MC = MonotoneComparator;

    ASSERT_TRUE(MC(BFDirection::Forward)({0, 0}, {1, 1}));
    ASSERT_FALSE(MC(BFDirection::Forward)({1, 1}, {0, 0}));

    ASSERT_FALSE(MC(BFDirection::Forward)({0, 0}, {0, 0}));
    ASSERT_FALSE(MC(BFDirection::Backward)({0, 0}, {0, 0}));
}

TEST(LineTest, IncludesPoint) {
    // Includes direction:
    ASSERT_TRUE(Line({2, 3}, {1, 0}).includesPoint({3, 3}));

    // Includes origin:
    ASSERT_TRUE(Line({2, 3}, {1, 0}).includesPoint({2, 3}));

    // Includes point some distance along the line:
    ASSERT_TRUE(Line({2, 3}, {1, 0}).includesPoint({7, 3}));

    // Does not include point outside line:
    ASSERT_FALSE(Line({2, 3}, {1, 0}).includesPoint({0, 1}));
}

TEST(LineTest, GetY) {
    // Diagonal:
    ASSERT_DOUBLE_EQ(Line({2, 1}, {1, 1}).getY(4), 3);

    // Vertical:
    ASSERT_TRUE(std::isnan(Line({2, 1}, {0, 1}).getY(4)));

    // Horizontal:
    ASSERT_DOUBLE_EQ(Line({2, 1}, {1, 0}).getY(4), 1);
}

TEST(LineTest, GetX) {
    // Diagonal:
    ASSERT_DOUBLE_EQ(Line({2, 1}, {1, 1}).getX(3), 4);

    // Vertical:
    ASSERT_DOUBLE_EQ(Line({2, 1}, {0, 1}).getX(4), 2);

    // Horizontal:
    ASSERT_TRUE(std::isnan(Line({2, 1}, {1, 0}).getX(4)));
}

TEST(LineTest, Intersect) {
    Point p = intersect(Line({0, 0}, {1, 1}), Line({0, 2}, {1, -1}));
    ASSERT_DOUBLE_EQ(p.x, 1);
    ASSERT_DOUBLE_EQ(p.y, 1);

    p = intersect(Line({0, 0}, {1, 1}), Line({0, 1}, {1, 0}));
    ASSERT_DOUBLE_EQ(p.x, 1);
    ASSERT_DOUBLE_EQ(p.y, 1);
}