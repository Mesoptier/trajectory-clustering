#include "simplification/agarwal.h"

#include <algorithm>

namespace {
    /**
     * \brief Get the area of the bounding box of a curve to use as upper bound
     * for the threshold.
     * \param c The curve.
     * \return The area of the axis-aligned bounding box fitted to c.
     */
    distance_t bb_area(const Curve& c) {
        const auto& extreme = c.getExtremePoints();
        auto a = extreme.max_x - extreme.min_x;
        auto b = extreme.max_y - extreme.min_y;
        return a * b;
    }

    /**
     * \brief Run a naive version that examines each point one-by-one; might be
     * slower than simplify_exponential due to a lot of distance computations.
     * \param in The original curve.
     * \param delta The distance threshold.
     * \param less_than The decision version of the distance function between
     * two curves, returns true if d(c1, c2) < dist.
     * \return The simplified curve.
     */
    Curve simplify_naive(const Curve& in, distance_t delta, const
            std::function<bool(const Curve&, const Curve&, distance_t)>&
            less_than) {
        Curve simplified_curve({in.front()});
        Curve prefix_curve({in.front()});

        for (PointID id = 1; id < in.size() - 1; ++id) {
            const auto& point = in[id];
            prefix_curve.push_back(point);
            auto segment = Curve({prefix_curve.front(), prefix_curve.back()});
            bool fits = prefix_curve.size() > 2 ?
                less_than(segment, prefix_curve, delta) : true;
            if (!fits) {
                simplified_curve.push_back(point);
                prefix_curve = Curve({point});
            }
        }
        simplified_curve.push_back(in.back());

        return simplified_curve;
    }

    /**
     * \brief Run exponential search followed by binary search on the
     * appropriate segment; might be slower than simplify_naive on small curves
     * due to a lot of curve point management.
     * \param in The original curve.
     * \param delta The distance threshold.
     * \param less_than The decision version of the distance function between
     * two curves, returns true if d(c1, c2) < dist.
     * \return The simplified curve.
     */
    Curve simplify_exponential(const Curve& in, distance_t delta, const
            std::function<bool(const Curve&, const Curve&, distance_t)>&
            less_than) {
        Curve simplified_curve({in.front()});
        Curve prefix_curve({in.front()});
        Curve next_curve({in.front()});

        PointID start = 0, exp = 1;
        while (start + exp < in.size()) {
            PointID next = std::min(start + 2 * exp, in.size());
            for (PointID i = exp; i < next; ++i)
                next_curve.push_back(in[i]);
            if (less_than(Curve({next_curve.front(), next_curve.back()}),
                          next_curve, delta)) {
                prefix_curve = next_curve;
                exp *= 2;
            }
            else {
                next_curve = prefix_curve;
                PointID max = next - 1, min = exp;
                while (max > min) {
                    auto mid = min + (max - min) / 2;
                    for (PointID i = min; i <= mid; ++i)
                        next_curve.push_back(in[i]);
                    if (less_than(Curve({next_curve.front(), next_curve.back()}),
                            next_curve, delta)) {
                        prefix_curve = next_curve;
                        min = mid + 1;
                    }
                    else {
                        next_curve = prefix_curve;
                        max = mid - 1;
                    }
                }
                simplified_curve.push_back(prefix_curve.back());
                prefix_curve = Curve({simplified_curve.back()});
                next_curve = prefix_curve;
                start += prefix_curve.size() - 1;
                exp = 1;
            }
        }

        simplified_curve.push_back(in.back());
        return simplified_curve;
    }
}

Curve simplification::greedy::simplify(const Curve& in, const PointID& ell,
        const std::function<bool(const Curve&, const Curve&, distance_t)>&
        less_than) {
    static constexpr distance_t epsilon = 1e-8;
    assert(in.size() > 1);
    assert(ell > 0);
    if (in.size() - 1 <= ell)
        return in;

    distance_t min = 0.0;
    distance_t max = bb_area(in);

    Curve simplified;
    while (max - min > epsilon) {
        auto split = (max + min) / 2.0;
        simplified = simplify(in, split, less_than);
        if (simplified.size() <= ell)
            max = split;
        else
            min = split;
    }

    return simplify(in, max, less_than);
}

Curve simplification::greedy::simplify(const Curve& in, distance_t delta, const
        std::function<bool(const Curve&, const Curve&, distance_t)>& less_than) {
    if (in.size() < 100)
        return simplify_naive(in, delta, less_than);
    return simplify_exponential(in, delta, less_than);
}
