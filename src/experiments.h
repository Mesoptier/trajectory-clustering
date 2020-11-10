#ifndef EXPERIMENTS_H
#define EXPERIMENTS_H

#include <cstddef>
#include <string>

namespace experiments {
    // Synthetic curves.
    void synthetic_curve_experiment();

    // Initial clustering (with pigeons).
    void initial_clustering_experiment();

    // Movebank.
    void center_update_experiment_movebank(std::string const& directory,
        std::size_t k, std::size_t l, bool remove_stops = false);

    // Characters.
    void center_update_experiment_characters(std::string const& directory,
        std::size_t n, std::size_t k, std::size_t l);

    void curve_complexity_experiment_characters();

    // Pigeons.
    void center_update_experiment_pigeons(std::string const& directory,
        std::size_t k, std::size_t l);

    void curve_complexity_experiment_pigeons();

    void find_wedge_params_pigeons();

    void test_windowed_convergence(unsigned step, unsigned w_size);
}
#endif
