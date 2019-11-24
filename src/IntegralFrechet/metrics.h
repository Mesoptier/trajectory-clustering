#ifndef CODE_METRICS_H
#define CODE_METRICS_H

enum class ImageMetric {
    L1,
    L2,
    L2_Squared,
};

enum class ParamMetric {
    L1,
    LInfinity_NoShortcuts,
};

#endif //CODE_METRICS_H
