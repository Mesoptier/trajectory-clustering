#include <geom.h>
#include "gtest/gtest.h"

TEST(PointTest, Perp) {
    // Collinear:
    ASSERT_DOUBLE_EQ(perp({1, 1}, {5, 5}), 0.0);

    // Second vector to the left of first:
    ASSERT_GT(perp({1, 1}, {0, 2}), 0);

    // Second vector to the right of first:
    ASSERT_LT(perp({1, 1}, {2, 0}), 0);
}

TEST(MonotoneComparatorTest, GetDirection) {
    using MC = MonotoneComparator;

    // Lower first:
    ASSERT_EQ(MC::getDirection({0, 0}, {1, 1}), MC::LowerFirst);
    // Higher first:
    ASSERT_EQ(MC::getDirection({1, 1}, {0, 0}), MC::HigherFirst);
    // Not monotone:
    ASSERT_ANY_THROW(MC::getDirection({1, 0}, {0, 1}));
    // Equal:
    ASSERT_ANY_THROW(MC::getDirection({0, 0}, {0, 0}));
}

TEST(MonotoneComparatorTest, Compare) {
    using MC = MonotoneComparator;

    ASSERT_TRUE(MC(MC::LowerFirst)({0, 0}, {1, 1}));
    ASSERT_FALSE(MC(MC::LowerFirst)({1, 1}, {0, 0}));

    // TODO: This is a perhaps unexpected case: equal points are both lower and higher.
    //  Especially odd since getDirection throws an error for equal points.
    ASSERT_TRUE(MC(MC::LowerFirst)({0, 0}, {0, 0}));
    ASSERT_TRUE(MC(MC::HigherFirst)({0, 0}, {0, 0}));
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