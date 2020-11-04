#include "simplification_experiment.h"

#include <array>
#include <chrono>
#include <cmath>
#include <exception>
#include <filesystem>
#include <iostream>
#include <numeric>
#include <sstream>
#include <type_traits>
#include <utility>

#include "distance_functions.h"
#include "simplification/agarwal.h"
#include "simplification/imaiiri.h"
#include "utils/io.h"

namespace {
    template<typename Duration = std::chrono::nanoseconds,
        typename F, typename ... Args>
    std::pair<std::invoke_result_t<F, Args...>, Duration>
            time_and_save(F&& f, Args&&... args) {
        auto const start = std::chrono::high_resolution_clock::now();
        auto res = std::forward<F>(f)(std::forward<Args>(args)...);
        auto const end = std::chrono::high_resolution_clock::now();
        return {res, std::chrono::duration_cast<Duration>(end - start)};
    }
}

using cvsz = std::vector<Curve>::size_type;

std::vector<Curve> experiments::sample(std::vector<Curve> const& curves,
        unsigned period, unsigned length) {
    std::vector<Curve> ret;
    for (cvsz i = 0; i < curves.size(); i += period)
        ret.emplace_back(curves[i].naive_l_simplification(length));
    return ret;
}

void experiments::evaluate(std::vector<Curve> const& curves,std::size_t ell,
        std::string const& prefix, unsigned scale) {
    using ms = std::chrono::milliseconds;

    std::array lt {df::dtw_lt, df::frechet_lt, df::integral_frechet_lt,
        df::integral_frechet_fast_lt};
    std::array dist {df::dtw, df::frechet, df::integral_frechet,
        df::integral_frechet_fast};
    std::array lt_names {"DTW    ", "Frechet", "CDTW   ", "F CDTW "};
    std::array lt_brief {"dtw", "fr", "cdtw", "fcdtw"};

    auto simpl1 = static_cast<Curve(*)(Curve const&, std::size_t,
        std::function<bool(Curve const&, Curve const&, distance_t)> const&)>
        (simplification::greedy::simplify);
    auto simpl2 = static_cast<Curve(*)(Curve const&, std::size_t,
        std::function<bool(Curve const&, Curve const&, distance_t)> const&)>
        (simplification::imai_iri::simplify);
    auto simpl_alt = [](Curve const& in, std::size_t l,
            std::function<distance_t(Curve const&, Curve const&)> const& d) {
        return simplification::imai_iri::simplify(in, l, d, false).second;
    };
    std::array algs {simpl1, simpl2};
    std::array alg_names {"Greedy  ", "Imai-Iri"};
    std::array alg_brief {"gr", "ii"};

    std::filesystem::path simpl_dir("simpl");
    simpl_dir /= prefix;
    try {
        std::filesystem::create_directories(simpl_dir);
    } catch (std::exception const& e) {
        std::cerr << e.what();
        return;
    }

    std::ostringstream stats;

    std::cout << std::string(20, ' ')
        << "time      mean       min       max    stddev\n";

    for (std::size_t i = 0; i < lt_names.size(); ++i) {
        auto const& lt_name = lt_names[i];
        auto const& ltf =  lt[i];
        auto const& distf = dist[i];

        for (std::size_t j = 0; j < algs.size(); ++j) {
            bool special = j == 1 && i != 1;
            auto const& alg = algs[j];
            auto const& alg_name = alg_names[j];
            std::vector<distance_t> cdtw_costs;
            std::chrono::nanoseconds time(0);

            for (auto const& c: curves) {
                auto r = special ? time_and_save(simpl_alt, c, ell, distf)
                                 : time_and_save(alg, c, ell, ltf);
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
                << std::setw(7) << std::chrono::duration_cast<ms>(time).count()
                << " " << (scale >= 7 ? std::scientific : std::fixed)
                << std::setprecision(scale < 7 ? 8 - static_cast<int>(scale) : 3)
                << std::setw(9) << mean_cost << " " << std::setw(9)
                << *extr.first << " " << std::setw(9) << *extr.second << " "
                << std::setw(9) << std::sqrt(var_est) << std::endl;

            if (j == 1 && (i == 2 || i == 3))
                continue;

            auto const worst_index = static_cast<std::size_t>(
                std::distance(cdtw_costs.begin(), extr.second));
            auto const& wc = curves[worst_index];
            auto r1 = special ? simpl_alt(wc, ell, distf) : alg(wc, ell, ltf);
            auto r2 = simpl_alt(wc, ell, df::integral_frechet);

            auto r1n = std::string(alg_brief[j]) + "_" + lt_brief[i];

            io::export_points(simpl_dir / ("crv" + std::to_string(worst_index)),
                wc.get_points());
            io::export_points(simpl_dir / (r1n + std::to_string(worst_index)),
                r1.get_points());
            io::export_points(simpl_dir / ("best" + std::to_string(worst_index)),
                r2.get_points());

            stats << "Cost for curve " << worst_index << ":\n"
                << "                CDTW cost Frechet cost\n" << std::fixed
                << std::setprecision(4) << std::setw(15)
                << alg_name << " + " << lt_name << " "
                << std::setw(9) << df::integral_frechet(wc, r1) << " "
                << std::setw(12) << df::frechet(wc, r1) << "\n"
                << "Imai-Iri + CDTW " << std::setw(9)
                << df::integral_frechet(wc, r2) << " " << std::setw(12)
                << df::frechet(wc, r2) << "\n" << std::endl;
        }
    }
    std::cout << "\n" << stats.str();
}
