#pragma once

#include <random>

class Random
{
public:
	Random();
	Random(unsigned int seed);

	unsigned int getSeed() const;

	int getUniformInt(int first, int last);
	bool throwCoin();

private:
	unsigned int seed;
	std::default_random_engine generator;
};
