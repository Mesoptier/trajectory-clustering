#include "clustering/pam.h"

#include <algorithm>
#include <limits>

#ifndef NDEBUG
#include <iostream>
#endif

namespace {
    template<class T>
    bool contains(const std::vector<T>& v, const T& x) {
        return v.end() != std::find(v.begin(), v.end(), x);
    }

    void build(std::size_t n, std::size_t k,
            const DistanceMatrix<distance_t>& d, distance_t& td,
            std::vector<std::size_t>& medoids) {
        // 1.
        td = std::numeric_limits<distance_t>::infinity();
        medoids.push_back(std::numeric_limits<std::size_t>::max());

        // 2.
        for (std::size_t xj = 0; xj < n; ++xj) {
            // 3.
            distance_t td_j = 0;
            // 4.
            for (std::size_t xo = 0; xo < n; ++xo) {
                if (xo != xj)
                    td_j += d.at(xo, xj);
            }
            // 5.
            if (td_j < td) {
                td = td_j;
                medoids[0] = xj;
            }
        }

        // Cache distance to nearest medoid
        std::vector<distance_t> d_nearest(n);
        for (std::size_t x = 0; x < n; ++x) {
            d_nearest[x] = d.at(x, medoids[0]);
        }

        // 6.
        for (std::size_t i = 1; i < k; ++i) {
            // 7.
            auto delta_td_best = std::numeric_limits<distance_t>::infinity();
            auto x_best = std::numeric_limits<std::size_t>::max();
            // 8.
            for (std::size_t xj = 0; xj < n; ++xj) {
                if (contains(medoids, xj))
                    continue;
                // 9.
                distance_t delta_td = 0;
                // 10.
                for (std::size_t xo = 0; xo < n; ++xo) {
                    if (xo == xj || contains(medoids, xo))
                        continue;
                    // 11.
                    distance_t delta = d.at(xo, xj) - d_nearest.at(xo);
                    // 12.
                    if (delta < 0) {
                        delta_td += delta;
                    }
                }
                // 13.
                if (delta_td < delta_td_best) {
                    delta_td_best = delta_td;
                    x_best = xj;
                }
            }
            // 14.
            td += delta_td_best;
            medoids.push_back(x_best);

            // Update cache
            for (std::size_t x = 0; x < n; ++x) {
                d_nearest[x] = std::min(d_nearest[x], d.at(x, x_best));
            }
        }
    }

    void swap(std::size_t n, const DistanceMatrix<distance_t>& d,
            distance_t& td, std::vector<std::size_t>& medoids) {
        std::vector<std::size_t> nearest(n);
        std::vector<distance_t> d_nearest(n,
            std::numeric_limits<distance_t>::infinity());
        std::vector<std::size_t> second_nearest(n);
        std::vector<distance_t> d_second_nearest(n,
            std::numeric_limits<distance_t>::infinity());

        for (std::size_t x = 0; x < n; ++x) {
            for (std::size_t mi : medoids) {
                distance_t dist = d.at(x, mi);
                if (dist < d_nearest[x]) {
                    second_nearest[x] = nearest[x];
                    d_second_nearest[x] = d_nearest[x];
                    nearest[x] = mi;
                    d_nearest[x] = dist;
                } else if (dist < d_second_nearest[x]) {
                    second_nearest[x] = mi;
                    d_second_nearest[x] = dist;
                }
            }
        }

        // 1.
        while (true) {
            // 2.
            auto delta_td_best = std::numeric_limits<distance_t>::infinity();
            auto m_best = std::numeric_limits<std::size_t>::max();
            auto x_best = std::numeric_limits<std::size_t>::max();
            // 3.
            for (std::size_t mi : medoids) {
                // 4.
                for (std::size_t xj = 0; xj < n; ++xj) {
                    if (contains(medoids, xj))
                        continue;
                    // 5.
                    distance_t delta_td = 0;
                    // 6.
                    for (std::size_t xo = 0; xo < n; ++xo) {
                        if (xo != mi && contains(medoids, xo))
                            continue;
                        distance_t delta = (mi == nearest.at(xo))
                            ? std::min(d.at(xo, xj), d_second_nearest.at(xo))
                                - d_nearest.at(xo)
                            : std::min(d.at(xo, xj) - d_nearest.at(xo), 0.0);
                        delta_td += delta;
                    }
                    // 7.
                    if (delta_td < delta_td_best) {
                        delta_td_best = delta_td;
                        m_best = mi;
                        x_best = xj;
                    }
                }
            }
            // 8.
            if (delta_td_best >= 0) {
                break;
            }
            // 9.
            std::replace(medoids.begin(), medoids.end(), m_best, x_best);

            // Update caches
            // TODO: Be more clever about this, instead of just brute-forcing it
            for (std::size_t x = 0; x < n; ++x) {
                d_nearest[x] = std::numeric_limits<distance_t>::infinity();
                d_second_nearest[x] = std::numeric_limits<distance_t>::infinity();

                for (std::size_t mi : medoids) {
                    distance_t dist = d.at(x, mi);
                    if (dist < d_nearest[x]) {
                        second_nearest[x] = nearest[x];
                        d_second_nearest[x] = d_nearest[x];
                        nearest[x] = mi;
                        d_nearest[x] = dist;
                    } else if (dist < d_second_nearest[x]) {
                        second_nearest[x] = mi;
                        d_second_nearest[x] = dist;
                    }
                }
            }

            // 10.
            td += delta_td_best;
        }
    }
}

std::vector<std::size_t> clustering::pam::compute(std::size_t n, std::size_t k,
        const DistanceMatrix<distance_t>& d) {
    distance_t td;
    std::vector<std::size_t> medoids;

    build(n, k, d, td, medoids);
    swap(n, d, td, medoids);

    std::sort(medoids.begin(), medoids.end());

    #ifndef NDEBUG
    std::cerr << "td: " << td << '\n';
    std::cerr << "medoids:";
    for (const auto& m: medoids)
        std::cerr << ' ' << m;
    std::cerr << std::endl;
    #endif

    return medoids;
}
