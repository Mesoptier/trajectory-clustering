#include "random.h"

#include <chrono>

Random::Random()
	: Random(std::chrono::system_clock::now().time_since_epoch().count())
{

}

Random::Random(unsigned int seed)
	: seed(seed)
{
	generator.seed(seed);
}

unsigned int Random::getSeed() const
{
	return seed;
}

int Random::getUniformInt(int first, int last)
{
	std::uniform_int_distribution<int> distribution(first, last);
	return distribution(generator);
}

bool Random::throwCoin()
{
	std::bernoulli_distribution distribution(0.5);
	return distribution(generator);
}
