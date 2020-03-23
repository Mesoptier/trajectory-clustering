#pragma once

#include "Curve.h"
#include "Frechet/frechet_light.h"



Curve simplify(Curve const& curve, int l, distance_t(*dist_func)(Curve, Curve));
