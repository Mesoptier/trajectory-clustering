#pragma once

#include "Polynomial.h"
#include "Interval.h"

enum CONST_TYPE { X_RIGHT, X_LEFT, Y_RANGE };

struct Constraint {
    Interval y_range;
    Polynomial<1> poly;
    CONST_TYPE type;
};

bool feasible_triple() {
    return true;
};