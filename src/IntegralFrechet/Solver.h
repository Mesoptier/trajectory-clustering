#ifndef CODE_SOLVER_H
#define CODE_SOLVER_H

#include <vector>
#include "Cell.h"
#include "../Curve.h"
#include "../expressionml.h"

template<class V>
class Solver
{
    Curve<double> curve1;
    Curve<double> curve2;

    // Curve complexity
    unsigned int n1;
    unsigned int n2;

    ImageMetric imageMetric;
    ParamMetric paramMetric;

    // Grid of cells
    std::vector<Cell<V>> cells;

public:
    Solver(const Curve<V>& curve1, const Curve<V>& curve2, double h, ImageMetric imageMetric, ParamMetric paramMetric);

    V getDistance() const;
    arma::Mat<V> getMatching() const;
    arma::Mat<V> getBoundaryCosts() const;

    void writeExpressionML(ExpressionML::Writer& writer) const {
        writer.openFunction("Association");

        // Image/parameter space metrics
        writer.openRule("ImageMetric");
        switch (imageMetric) {
            case ImageMetric::L1:
                writer.writeString("L1");
                break;
            case ImageMetric::L2:
                writer.writeString("L2");
                break;
            case ImageMetric::L2_Squared:
                writer.writeString("L2_Squared");
                break;
        }
        writer.closeRule();

        // Cells
        writer.openRule("Cells");
        writer.openFunction("List");
        for (int i1 = 0; i1 < n1; ++i1) {
            writer.openFunction("List");
            for (int i2 = 0; i2 < n2; ++i2) {
                cells[i1 * n2 + i2].writeExpressionML(writer);
            }
            writer.closeFunction();
        }
        writer.closeFunction();
        writer.closeRule();

        writer.closeFunction();
    }
};

#endif //CODE_SOLVER_H
