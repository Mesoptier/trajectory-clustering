#include "Frechet/frechet_matching.h"

#include "Frechet/frechet_light.h"
#include "utils/defs.h"

Points frechet::calcMatching(Curve const& curve1, Curve const& curve2) {
    FrechetLight frechet_light;

    auto dist = frechet_light.calcDistance(curve1, curve2);
    frechet_light.lessThan(dist*(1+10e-4), curve1, curve2);

    Points matching;
    frechet_light.computeCertificate();
    auto cert = frechet_light.getCertificate();
    if (!cert.isYes()) {
        ERROR("The certificate should be positive!");
    }

    matching.push_back(curve2.front());
    for (auto const& cpos: cert.getTraversal()) {
        auto cpoint1 = cpos[0];
        auto cpoint2 = cpos[1];
        while (matching.size() <= cpoint1.getPoint())
            matching.push_back(curve2.interpolate_at(cpoint2));
    }
    matching.pop_back();
    matching.push_back(curve2.back());
    return matching;
}
