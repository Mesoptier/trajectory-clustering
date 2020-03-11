#include <cmath>
#include <algorithm>
#include <iostream>
#include <map>
#include "pam.h"


namespace {
    template<class T>
    bool contains(const std::vector<T>& v, const T& x) {
        return v.end() != std::find(v.begin(), v.end(), x);
    }

    void build(size_t n, size_t k, const SymmetricMatrix& d, double& td, std::vector<size_t>& medoids) {
        // 1.
        td = INFINITY;
        medoids.push_back(std::numeric_limits<size_t>::max());

        // 2.
        for (size_t xj = 0; xj < n; ++xj) {
            // 3.
            double td_j = 0;
            // 4.
            for (size_t xo = 0; xo < n; ++xo) {
                if (xo == xj) continue;
                td_j += d.at(xo, xj);
            }
            // 5.
            if (td_j < td) {
                td = td_j;
                medoids[0] = xj;
            }
        }

        // Cache distance to nearest medoid
        std::vector<double> d_nearest(n);
        for (size_t x = 0; x < n; ++x) {
            d_nearest[x] = d.at(medoids[0], x);
        }

        // 6.
        for (size_t i = 1; i < k; ++i) {
            // 7.
            double delta_td_best = INFINITY;
            size_t x_best = std::numeric_limits<size_t>::max();
            // 8.
            for (size_t xj = 0; xj < n; ++xj) {
                if (contains(medoids, xj)) continue;
                // 9.
                double delta_td = 0;
                // 10.
                for (size_t xo = 0; xo < n; ++xo) {
                    if (xo == xj || contains(medoids, xo)) continue;
                    // 11.
                    double delta = d.at(xo, xj) - d_nearest.at(xo);
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
            for (size_t x = 0; x < n; ++x) {
                d_nearest[x] = std::min(d_nearest[x], d.at(x_best, x));
            }
        }
    }

    void swap(size_t n, size_t k, const SymmetricMatrix& d, double& td, std::vector<size_t>& medoids) {
        std::vector<size_t> nearest(n);
        std::vector<double> d_nearest(n, INFINITY);
        std::vector<size_t> second_nearest(n);
        std::vector<double> d_second_nearest(n, INFINITY);

        for (size_t x = 0; x < n; ++x) {
            for (size_t mi : medoids) {
                double dist = d.at(x, mi);
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
            double delta_td_best = INFINITY;
            size_t m_best = std::numeric_limits<size_t>::max();
            size_t x_best = std::numeric_limits<size_t>::max();
            // 3.
            for (size_t mi : medoids) {
                // 4.
                for (size_t xj = 0; xj < n; ++xj) {
                    if (contains(medoids, xj)) continue;
                    // 5.
                    double delta_td = 0;
                    // 6.
                    for (size_t xo = 0; xo < n; ++xo) {
                        if (xo != mi && contains(medoids, xo)) continue;
                        double delta = (mi == nearest.at(xo))
                            ? std::min(d.at(xo, xj), d_second_nearest.at(xo)) - d_nearest.at(xo)
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
            for (size_t x = 0; x < n; ++x) {
                d_nearest[x] = INFINITY;
                d_second_nearest[x] = INFINITY;

                for (size_t mi : medoids) {
                    double dist = d.at(x, mi);
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

        //
        // EXPORT MEDOIDS
        // TODO: Return medoids instead
        //
        std::map<size_t, size_t> medoid_indices;
        for (size_t i = 0; i < k; ++i) {
            medoid_indices[medoids[i]] = i;
        }

        SymmetricMatrix clusters(n);
        for (size_t xi = 0; xi < n; ++xi) {
            for (size_t xj = 0; xj < n; ++xj) {
                clusters.at(xi, xj) = (nearest[xi] == nearest[xj]) ? medoid_indices.at(nearest[xi]) : 0;
            }
        }
        std::ofstream file("data/out/clusters.mtx");
        clusters.write(file);
        file.close();
    }
}

namespace clustering::pam {

    void compute(size_t n, size_t k, const SymmetricMatrix& d) {
        double td;
        std::vector<size_t> medoids;

        build(n, k, d, td, medoids);
        swap(n, k, d, td, medoids);

        std::sort(medoids.begin(), medoids.end());

        std::cout << "td: " << td << '\n';
        std::cout << "medoids: ";
        for (auto m : medoids) {
            std::cout << m << ' ';
        }
        std::cout << std::endl;
    }
}

