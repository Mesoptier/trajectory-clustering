#include "gtest/gtest.h"
#include "../src/cdtw/Polynomial.h"

::testing::AssertionResult IsApproxEqual(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) {
        return ::testing::AssertionFailure() << "different sizes";
    }

    for (size_t i = 0; i < a.size(); ++i) {
        if (!approx_equal(a[i], b[i])) {
            return ::testing::AssertionFailure() << a[i] << " !~= " << b[i];
        }
    }

    return ::testing::AssertionSuccess();
}

TEST(PolynomialTest, FindRootsCubic) {
    {
        // 3 real roots
        Polynomial<3> f({0, -1, 0, 1 });
        const auto actual = find_roots(f);
        const std::vector<double> expected{1, 0, -1};
        ASSERT_TRUE(IsApproxEqual(actual, expected));
    }
    {
        // 1 real triple root
        Polynomial<3> f({0, 0, 0, 1 });
        const auto actual = find_roots(f);
        const std::vector<double> expected{0};
        ASSERT_TRUE(IsApproxEqual(actual, expected));
    }
    {
        // 1 real root + 1 real double root
        Polynomial<3> f({0, 0, 1, 1 });
        const auto actual = find_roots(f);
        const std::vector<double> expected{-1, 0};
        ASSERT_TRUE(IsApproxEqual(actual, expected));
    }
}