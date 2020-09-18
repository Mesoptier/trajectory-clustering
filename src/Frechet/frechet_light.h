#ifndef FRECHET_LIGHT_H
#define FRECHET_LIGHT_H

#include <array>
#include <vector>
#include "Curve.h"
#include "Frechet/certificate.h"
#include "Frechet/filter.h"
#include "Frechet/frechet_abstract.h"
#include "Frechet/frechet_light_types.h"
#include "geom.h"
#include "utils/defs.h"
#include "utils/id.h"

class FrechetLight final : public FrechetAbstract {
    using CurvePair = std::array<Curve const*, 2>;

public:
    static constexpr distance_t eps = 1e-10;

    void buildFreespaceDiagram(distance_t distance, Curve const& curve1,
        Curve const& curve2);
    bool lessThan(distance_t distance, Curve const& curve1,
        Curve const& curve2) override;
    bool lessThanWithFilters(distance_t distance, Curve const& curve1,
        Curve const& curve2);
    distance_t calcDistance(Curve const& curve1, Curve const& curve2);
    void clear();

    CurvePair getCurvePair() const;
    Certificate& computeCertificate() override;
    Certificate const& getCertificate() const {
        return cert;
    }

    void setPruningLevel(int pruning_level) override;
    void setRules(std::array<bool, 5> const& enable) override;

    std::size_t getNumberOfBoxes() const;

    std::size_t non_filtered = 0;

private:
    CurvePair curve_pair;

    distance_t distance;
    distance_t dist_sqr;

    std::vector<CIntervals> reachable_intervals_vec;
    frechet::QSimpleIntervals qsimple_intervals;
    std::size_t num_boxes;

    // 0 = no pruning ... 6 = full pruning
    int pruning_level = 6;
    // ... and additionally bools to enable/disable rules
    bool enable_box_shrinking = true;
    bool enable_empty_outputs = true;
    bool enable_propagation1 = true;
    bool enable_propagation2 = true;
    bool enable_boundary_rule = true;

#ifdef VIS
    CIntervals unknown_intervals;
    CIntervals connections;
    CIntervals free_non_reachable;
    CIntervals reachable_intervals;

    struct Cell {
        PointID i, j;
        Cell(PointID i, PointID j) : i(i), j(j) {}
    };

    std::vector<Cell> cells;
    std::vector<Cell> const& getCells() const {
        return cells;
    }
#endif

    Certificate cert;

    CInterval getInterval(Point const& point, Curve const& curve,
        PointID i) const;
    CInterval getInterval(Point const& point, Curve const& curve,
        PointID i, CInterval* ) const;
    void merge(CIntervals& v, CInterval const& i) const;

    frechet::Outputs createFinalOutputs();
    frechet::Inputs computeInitialInputs();
    // XXX: consistency of arguments in following functions!
    distance_t getDistToPointSqr(Curve const& curve, Point const& point) const;
    bool isClose(Point const& point, Curve const& curve) const;
    CPoint getLastReachablePoint(Point const& point, Curve const& curve) const;
    bool isTopRightReachable(frechet::Outputs const& outputs) const;
    void computeOutputs(frechet::Box const& initial_box,
        frechet::Inputs const& initial_inputs, frechet::Outputs& final_outputs);

    void getReachableIntervals(frechet::BoxData& data);

    // subfunctions of getReachableIntervals
    bool emptyInputsRule(frechet::BoxData& data);
    void boxShrinkingRule(frechet::BoxData& data);
    void handleCellCase(frechet::BoxData& data);
    void getQSimpleIntervals(frechet::BoxData& data);
    void calculateQSimple1(frechet::BoxData& data);
    void calculateQSimple2(frechet::BoxData& data);
    bool boundaryPruningRule(frechet::BoxData& data);
    void splitAndRecurse(frechet::BoxData& data);

    // intervals used in getReachableIntervals and subfunctions
    CInterval const empty;
    CInterval const* firstinterval1;
    CInterval const* firstinterval2;
    distance_t min1_frac, min2_frac;
    frechet::QSimpleInterval qsimple1, qsimple2;
    CInterval out1, out2;
    // TODO: can those be made members of out1, out2?
    bool out1_valid = false, out2_valid = false;

    // qsimple interval calculation functions
    frechet::QSimpleInterval getFreshQSimpleInterval(Point const& fixed_point,
        PointID min1, PointID max1, Curve const& curve) const;
    bool updateQSimpleInterval(frechet::QSimpleInterval& qsimple,
        Point const& fixed_point, PointID min1, PointID max1,
        Curve const& curve) const;
    void continueQSimpleSearch(frechet::QSimpleInterval& qsimple,
        Point const& fixed_point, PointID min1, PointID max1,
        Curve const& curve) const;

    bool isOnLowerRight(CPosition const& pt) const;
    bool isOnUpperLeft(CPosition const& pt) const;

    void initCertificate(frechet::Inputs const& initial_inputs);
    void certSetValues(CInterval& interval, CInterval const& parent,
        PointID point_id, CurveID curve_id);
    void certAddEmpty(CPoint begin, CPoint end, CPoint fixed_point,
        CurveID fixed_curve);

    // Those are empty function if VIS is not defined
    void visAddReachable(CInterval const& cinterval);
    void visAddUnknown(CPoint begin, CPoint end, CPoint fixed_point,
        CurveID fixed_curve);
    void visAddConnection(CPoint begin, CPoint end, CPoint fixed_point,
        CurveID fixed_curve);
    void visAddFreeNonReachable(CPoint begin, CPoint end, CPoint fixed_point,
        CurveID fixed_curve);
    void visAddCell(frechet::Box const& box);

    // Could also be done via getter member functions, but vis is a special
    // case in needing access to the internal structures.
    friend class FreespaceLightVis;
};
#endif
