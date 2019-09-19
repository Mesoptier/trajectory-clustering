#include "Curve.h"
#include "FastMarchCDTW.h"

int main() {
    Curve<double> curve1({{0, 0},
                    {1, 0}});

    Curve<double> curve2({{0, 1},
                          {1, 1}});

    std::vector<double> hs = {
        1,
        .1,
        .01,
        .001,
    };

    for (auto h : hs) {
        std::cout << h << '\t' << FastMarchCDTW::compute(curve1, curve2, h) << std::endl;
    }

    return 0;
}