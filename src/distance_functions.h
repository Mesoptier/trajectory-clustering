#pragma once

#ifndef DISTANCE_FUNCTIONS_H
#define DISTANCE_FUNCTIONS_H

#include "Curve.h"

/**
 * \brief Convenient definitions of all the distance metrics used.
 */
namespace df {
    /**
     * \brief Compute the dynamic time warping distance between two curves.
     * \param curve_1 First curve.
     * \param curve_2 Second curve.
     * \return The distance.
     */
    distance_t dtw(Curve const& curve_1, Curve const& curve_2);

    /**
     * \brief Compute the CDTW distance between two curves (slow approximation).
     * \param curve_1 First curve.
     * \param curve_2 Second curve.
     * \return The distance.
     */
    distance_t integral_frechet(Curve const& curve_1, Curve const& curve_2);

    /**
     * \brief Compute the CDTW distance between two curves (fast approximation).
     *
     * Should be fine for `normal' curves, but for weird curves it might
     * overestimate the actual distance by not considering some couplings.
     * \param curve_1 First curve.
     * \param curve_2 Second curve.
     * \return The distance.
     */
    distance_t integral_frechet_fast(Curve const& curve_1,
        Curve const& curve_2);

    /**
     * \brief Compute the Fréchet distance between two curves.
     * \param curve_1 First curve.
     * \param curve_2 Second curve.
     * \return The distance.
     */
    distance_t frechet(Curve const& curve_1, Curve const& curve_2);

    /**
     * \brief Compute the average Fréchet distance between two curves
     * (integral_frechet) with curve-length adjustments.
     * \param curve_1 First curve.
     * \param curve_2 Second curve.
     * \return The distance.
     */
    distance_t average_frechet(Curve const& curve_1, Curve const& curve_2);


    /**
     * @brief Compute the 2d_l1_l1 cdtw distance vetween two curves 
     * 
     * @param curve_1 First curve.
     * @param curve_2 Second curve.
     * @return The distance. 
     */
    distance_t cdtw_2d_l1_l1(Curve const& curve_1, Curve const& curve_2);

    /**
     * @brief Compute the 2d_l1_l1 cdtw distance vetween two curves using the heurisitic algorithm 
     * 
     * @param curve_1 First curve.
     * @param curve_2 Second curve.
     * @return The distance. 
     */
    distance_t heur_cdtw_2d_l1_l1(Curve const& curve_1, Curve const& curve_2);

    /**
     * @brief Compute the 2d_l1_l1 cdtw distance vetween two curves using the heurisitic algorithm 
     * and write the resulting warping to the file warping_path_h.txt
     * 
     * @param curve_1 First curve.
     * @param curve_2 Second curve.
     * @return The distance.
     */
    distance_t heur_cdtw_2d_l1_l1_path(Curve const& curve_1, Curve const& curve_2);

    /**
     * \brief Compute the average Fréchet distance between two curves
     * (integral_frechet_fast) with curve-length adjustments.
     * \param curve_1 First curve.
     * \param curve_2 Second curve.
     * \return The distance.
     */
    distance_t average_frechet_fast(Curve const& curve_1, Curve const& curve_2);

    /**
     * \brief Check if the dynamic time warping distance between two curves is
     * below the threshold.
     * \param curve_1 First curve.
     * \param curve_2 Second curve.
     * \param delta The distance threshold.
     * \return True iff dtw(curve_1, curve_2) <= delta.
     */
    bool dtw_lt(Curve const& curve_1, Curve const& curve_2, distance_t delta);

    /**
     * \brief Check if the CDTW distance (slow approximation) between two curves
     * is below the threshold.
     * \param curve_1 First curve.
     * \param curve_2 Second curve.
     * \param delta The distance threshold.
     * \return True iff integral_frechet(curve_1, curve_2) <= delta.
     */
    bool integral_frechet_lt(Curve const& curve_1, Curve const& curve_2,
        distance_t delta);

    /**
     * \brief Check if the CDTW distance (fast approximation) between two curves
     * is below the threshold.
     * \param curve_1 First curve.
     * \param curve_2 Second curve.
     * \param delta The distance threshold.
     * \return True iff integral_frechet_fast(curve_1, curve_2) <= delta.
     */
    bool integral_frechet_fast_lt(Curve const& curve_1, Curve const& curve_2,
        distance_t delta);

    /**
     * \brief Check if the Fréchet distance between two curves is below the
     * threshold.
     * \param curve_1 First curve.
     * \param curve_2 Second curve.
     * \param delta The distance threshold.
     * \return True iff d_F(curve_1, curve_2) <= delta.
     */
    bool frechet_lt(Curve const& curve_1, Curve const& curve_2,
        distance_t delta);

    /**
     * \brief Check if the average Fréchet distance (slow) between two curves is
     * below the threshold.
     * \param curve_1 First curve.
     * \param curve_2 Second curve.
     * \param delta The distance threshold.
     * \return True iff average_frechet(curve_1, curve_2) <= delta.
     */
    bool average_frechet_lt(Curve const& curve_1, Curve const& curve_2,
        distance_t delta);

    /**
     * \brief Check if the average Fréchet distance (fast) between two curves is
     * below the threshold.
     * \param curve_1 First curve.
     * \param curve_2 Second curve.
     * \param delta The distance threshold.
     * \return True iff average_frechet_fast(curve_1, curve_2) <= delta.
     */
    bool average_frechet_fast_lt(Curve const& curve_1, Curve const& curve_2,
        distance_t delta);

}
#endif
