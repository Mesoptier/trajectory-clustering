#ifndef PLOT_CLUSTERING_H
#define PLOT_CLUSTERING_H

#include "basic_types.h"

namespace clustering {
    void plot_clustering(Clustering const& clustering,
        std::vector<Curve> const& curves, std::string const& filename);
}
#endif
