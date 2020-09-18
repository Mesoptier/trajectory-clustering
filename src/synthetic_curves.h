#ifndef SYNTHETIC_CURVES_H
#define SYNTHETIC_CURVES_H

#include <random>

#include "Curve.h"
#include "utils/io.h"

namespace synth {
    /**
     * \brief Generate displacement vectors for points.
     */
    class DisplacementGenerator {
        std::mt19937_64 rng;
        double lower;
        double upper;

    public:
        /**
         * \brief Set up the generator.
         * \param lower_lim Lower bound limit for displacement magnitude.
         * \param upper_lim Upper bound limit for displacement magnitude.
         * \param seed The random generator seed (if it needs to be fixed).
         */
        DisplacementGenerator(double lower_lim, double upper_lim,
            std::mt19937_64::result_type seed = std::random_device()());

        /**
         * \brief Generate a displacement vector in 2D. Any angle; magnitude
         * between lower_lim and upper_lim.
         * \return x and y coordinates of displacement vector.
         */
        Point generate_disp_vec();
    };

    /**
     * \brief Using the base curve, generate a certain number of randomly
     * displaced curves. All are translated; individual points are also
     * translated.
     * \param curve The true curve.
     * \param count The number of curves to generate.
     * \return The generated curves.
     */
    std::vector<Curve> generate_curves(Curve const& curve, std::size_t count);

    /**
     * \brief Generate 20 curves starting from base curve and save them in files
     * in a subdirectory.
     * \param base The true curve.
     */
    void write_curves(Curve const& base);
}
#endif
