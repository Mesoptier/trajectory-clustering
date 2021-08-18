// #include "gtest/gtest.h"
// #include "../src/Curve.h"
// #include "../src/cdtw/cdtw.h"
// #include "../src/cdtw/1d-l2squared-l1.h"

// using MyCDTW = CDTW<1, Norm::L2Squared, Norm::L1>;
// using MyFunctions = std::vector<std::vector<MyCDTW::Entry>>;

// void AssertFunctionsEqual(const MyFunctions& actual, const MyFunctions& expected) {
//     ASSERT_EQ(actual.size(), expected.size());
//     for (size_t i = 0; i < actual.size(); ++i) {
//         ASSERT_EQ(actual[i].size(), expected[i].size());
//         for (size_t j = 0; j < actual[i].size(); ++j) {
//             std::cout << "i: " << i << " j: " << j << '\n';
//             std::cout << "bottom:\n";
//             std::cout << "|   actual: " << actual[i][j].bottom << '\n';
//             std::cout << "| expected: " << expected[i][j].bottom << '\n';
//             std::cout << "left:\n";
//             std::cout << "|   actual: " << actual[i][j].left << '\n';
//             std::cout << "| expected: " << expected[i][j].left << '\n';
//             ASSERT_TRUE(approx_equal(actual[i][j].bottom, expected[i][j].bottom));
//             ASSERT_TRUE(approx_equal(actual[i][j].left, expected[i][j].left));
//         }
//     }
// }

// TEST(CDTW_L1_L2Squared_L1, Case1) {
//     const Curve curve1("curve1", {0, 1});
//     const Curve curve2("curve2", {1, 2});
//     const MyCDTW alg(curve1, curve2);
//     const MyFunctions& actual = alg.get_functions();
//     const MyFunctions expected{
//         {
//             {
//                 .bottom = PiecewisePolynomial<3>({
//                     {{0, 1}, Polynomial<3>({0, 1, -1, 1./3})},
//                 }),
//                 .left = PiecewisePolynomial<3>({
//                     {{0, 1}, Polynomial<3>({0, 1, 1, 1./3})},
//                 }),
//             },
//             {
//                 .bottom = PiecewisePolynomial<3>({
//                     {{0, 1}, Polynomial<3>({2 + 1./3, -2, 0, 1./3})},
//                 }),
//             },
//         },
//         {
//             {
//                 .left = PiecewisePolynomial<3>({
//                     {{0, 1}, Polynomial<3>({1./3, 0, 0, 1./3})},
//                 }),
//             },
//             {},
//         },
//     };

//     AssertFunctionsEqual(actual, expected);
// }

// TEST(CDTW_L1_L2Squared_L1, Case2) {
//     const Curve curve1("curve1", {0, 1});
//     const Curve curve2("curve2", {2, 1});
//     const MyCDTW alg(curve1, curve2);
//     const MyFunctions& actual = alg.get_functions();
//     const MyFunctions expected{
//         {
//             {
//                 .bottom = PiecewisePolynomial<3>({
//                     {{0, 1}, Polynomial<3>({0, 4, -2, 1./3})},
//                 }),
//                 .left = PiecewisePolynomial<3>({
//                     {{0, 1}, Polynomial<3>({0, 4, -2, 1./3})},
//                 }),
//             },
//             {
//                 .bottom = PiecewisePolynomial<3>({
//                     {{0, 1}, Polynomial<3>({2 + 1./3, 1, -1, 1./3})},
//                 }),
//             },
//         },
//         {
//             {
//                 .left = PiecewisePolynomial<3>({
//                     {{0, 1}, Polynomial<3>({2 + 1./3, 1, -1, 1./3})},
//                 }),
//             },
//             {},
//         },
//     };

//     AssertFunctionsEqual(actual, expected);
// }

// TEST(CDTW_L1_L2Squared_L1, Case3) {
//     const Curve curve1("curve1", {0, 1});
//     const Curve curve2("curve2", {0, 1});
//     const MyCDTW alg(curve1, curve2);
//     const MyFunctions& actual = alg.get_functions();
//     const MyFunctions expected{
//         {
//             {
//                 .bottom = PiecewisePolynomial<3>({
//                     {{0, 1}, Polynomial<3>({0, 0, 0, 1./3})},
//                 }),
//                 .left = PiecewisePolynomial<3>({
//                     {{0, 1}, Polynomial<3>({0, 0, 0, 1./3})},
//                 }),
//             },
//             {
//                 .bottom = PiecewisePolynomial<3>({
//                     {{0, 1}, Polynomial<3>({1./3, -1, 1, -1./3})},
//                 }),
//             },
//         },
//         {
//             {
//                 .left = PiecewisePolynomial<3>({
//                     {{0, 1}, Polynomial<3>({1./3, -1, 1, -1./3})},
//                 }),
//             },
//             {},
//         },
//     };

//     AssertFunctionsEqual(actual, expected);
// }

//TEST(CDTW_L1_L2Squared_L1, Case4) {
//    const Curve curve1("curve1", {0, 1});
//    const Curve curve2("curve2", {0, 1, 0});
//    const MyCDTW alg(curve1, curve2);
//    const MyFunctions& actual = alg.get_functions();
//    const MyFunctions expected{};
//
//    AssertFunctionsEqual(actual, expected);
//}
