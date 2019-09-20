#include "Curve.h"
#include "FastMarchCDTW.h"

int main() {
    Curve<double> curve1({{0, 0},
                          {0.5, 0.99},
                    {1, 0}});

    Curve<double> curve2({{0, 1},
                          {0.5, 1},
                          {1, 1}});

    FastMarchCDTW::compute(curve1, curve2, 0.1, true);

    return 0;
}