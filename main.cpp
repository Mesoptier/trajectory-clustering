#include <iostream>
#include <string>

#include "src/Curve.h"
#include "src/classification_experiment.h"
#include "src/experiments.h"
#include "src/simplification_experiment.h"
#include "src/utils/io.h"



int main() {

    experiments::characters_param_selection();
    experiments::aic();
    experiments::elbow_graph();
}
