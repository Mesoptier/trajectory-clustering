#pragma once

#include "../basic_types.h"
#include "../defs.h"

using Curves = std::vector<Curve>;

Curve cdba_characters(Curves const& curves, Cluster const& cluster);

Curve cdba_pigeons(Curves const& curves, Cluster const& cluster);

Curve cdba_storks(Curves const& curves, Cluster const& cluster);

Curve cdba_geo(Curves const& curves, Cluster const& cluster);

Curve dba_characters(Curves const& curves, Cluster const& cluster);

Curve dba_geo(Curves const& curves, Cluster const& cluster);

Curve dba_storks(Curves const& curves, Cluster const& cluster);

Curve dba_pigeons(Curves const& curves, Cluster const& cluster);

Curve wedge_characters(Curves const& curves, Cluster const& cluster);

Curve wedge_pigeons(Curves const& curves, Cluster const& cluster);

Curve wedge_storks(Curves const& curves, Cluster const& cluster);

Curve cdba_update(Curves const& curves, Cluster const& cluster);

Curve cdba_update(Curves const& curves, Cluster const& cluster, int res);

Curve wedge_update(Curves const& curves, Cluster const& cluster);

Curve regression_update(Curves const& curves, Cluster const& cluster);

Curve regression_update_3d(Curves const& curves, Cluster const& cluster);

Curve dba_update(Curves const& curves, Cluster const& cluster);

Curve fsa_update(Curves const& curves, Cluster const& cluster);

Curve redistribute_points_update(Curves const& curves, Cluster const& cluster);

Curve wedge_update_param_args(Curves const& curves, Cluster const& cluster, distance_t eps, int radius);
