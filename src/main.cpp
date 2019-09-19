#include "Curve.h"
#include "FastMarchCDTW.h"

int main() {
    Curve<double> curve1({{0, 0},
                    {1, 0},
                    {2, 0}});

    Curve<double> curve2({{0, 1},
                          {1, 1},
                          {2, 1}});

    FastMarchCDTW::compute(curve1, curve2);

    return 0;
}