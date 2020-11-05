#ifndef RANDOM_H
#define RANDOM_H

#include <random>
#include <type_traits>

/**
 * \brief Simple wrapper to generate random numbers and booleans.
 */
class Random
{
    std::mt19937_64 generator;
public:
    /**
     * \brief Get a random integer between two values.
     * \tparam T Integer type, one of short, int, long, long long, or the
     * unsigned versions.
     * \param first Lowest possible value.
     * \param last Highest possible value.
     * \return Uniformly chosen integer from [first, last].
     */
    template<typename T>
    std::enable_if_t<std::is_integral_v<T>, T> getUniformInt(T first, T last) {
        std::uniform_int_distribution<T> distribution(first, last);
        return distribution(generator);
    }

    /**
     * \brief Get the value from a Bernoulli distribution with p = 0.5.
     * \return True or False, with equal probability.
     */
    bool throwCoin() {
        std::bernoulli_distribution distribution(0.5);
        return distribution(generator);
    }
};
#endif
