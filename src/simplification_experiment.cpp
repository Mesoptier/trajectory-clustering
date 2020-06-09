#include "simplification_experiment.h"

#include <array>
#include <chrono>
#include <cmath>
#include <iostream>
#include <numeric>
#include <type_traits>
#include <utility>

#include "distance_functions.h"
#include "simplification/agarwal.h"
#include "simplification/imaiiri.h"

namespace {
    template<typename Duration = std::chrono::nanoseconds,
        typename F, typename ... Args>
    std::pair<std::invoke_result_t<F, Args...>, Duration>
            time_and_save(F&& f, Args&&... args) {
        const auto start = std::chrono::high_resolution_clock::now();
        auto res = std::forward<F>(f)(std::forward<Args>(args)...);
        const auto end = std::chrono::high_resolution_clock::now();
        return {res, std::chrono::duration_cast<Duration>(end - start)};
    }
}

using cvsz = std::vector<Curve>::size_type;

std::vector<Curve> experiments::sample(std::vector<Curve> const& curves,
        unsigned period) {
    std::vector<Curve> ret;
    for (cvsz i = 0; i < curves.size(); i += period)
        ret.emplace_back(curves[i]);
    return ret;
}

void experiments::evaluate(std::vector<Curve> const& curves,
        std::size_t ell) {
    using ms = std::chrono::milliseconds;

    std::array lt {df::dtw_lt, df::frechet_lt, df::integral_frechet_lt,
        df::integral_frechet_fast_lt};
    std::array lt_names {"DTW    ", "Frechet", "CDTW   ", "F CDTW "};

    auto simpl1 = static_cast<Curve(*)(Curve const&, PointID const&,
        std::function<bool(Curve const&, Curve const&, distance_t)> const&)>
        (simplification::greedy::simplify);
    auto simpl2 = static_cast<Curve(*)(Curve const&, PointID const&,
        std::function<bool(Curve const&, Curve const&, distance_t)> const&)>
        (simplification::imai_iri::simplify);
    std::array algs {simpl1, simpl2};
    std::array alg_names {"Greedy  ", "Imai-Iri"};

    std::cout << "                   time    mean     min     max  stddev\n";

    for (std::size_t i = 0; i < lt.size(); ++i) {
        auto const& ltf = lt[i];
        auto const& lt_name = lt_names[i];
        for (std::size_t j = 0; j < algs.size(); ++j) {
            auto const& alg = algs[j];
            auto const& alg_name = alg_names[j];
            std::vector<distance_t> cdtw_costs;
            std::chrono::nanoseconds time(0);

            if (j == 1 && i == 2) {
                std::cout << lt_name << " " << alg_name << "    ---     --- "
                    << "    ---     ---     ---\n";
                    continue;
            }

            for (auto const& c: curves) {
                auto r = time_and_save(alg, c, ell, ltf);
                cdtw_costs.emplace_back(df::integral_frechet(c, r.first));
                time += r.second;
            }

            auto const sz = curves.size();
            auto mean_cost =
                std::accumulate(cdtw_costs.begin(), cdtw_costs.end(), 0.0) / sz;
            auto extr = std::minmax_element(cdtw_costs.begin(),
                                            cdtw_costs.end());
            auto var_est = std::accumulate(cdtw_costs.begin(), cdtw_costs.end(),
                0.0, [&mean_cost, &sz](distance_t acc, distance_t const& v) {
                    return acc + (v - mean_cost) * (v - mean_cost) / (sz - 1);
            });
            std::cout << lt_name << " " << alg_name << " "
                << std::setw(6) << std::chrono::duration_cast<ms>(time).count()
                << " " << std::fixed << std::setw(7) << std::setprecision(3)
                << mean_cost << " " << std::setw(7) << *extr.first << " "
                << std::setw(7) << *extr.second << " " << std::setw(7)
                << std::sqrt(var_est) << std::endl;

        }
    }
}
