#include <iostream>
#include <armadillo>

using namespace arma;

int main() {
    mat A = randu<mat>(4, 5);

    std::cout << A << std::endl;

    return 0;
}