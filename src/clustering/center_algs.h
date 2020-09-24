#ifndef CENTER_ALGS_H
#define CENTER_ALGS_H

#include <functional>
#include <string>
#include <vector>

#include "basic_types.h"
#include "clustering/center_update.h"
#include "utils/defs.h"

namespace clustering {
    using Curves = std::vector<Curve>;

    /**
     * \brief Iterative center update methods (see paper).
     */
    enum class CenterAlg {
        fsa, dba, cdba, wedge, none
    };

    /**
     * \brief Get the name of the algorithm.
     * \param center_alg The algorithm used.
     * \return The name (FSA, DBA, CDBA, Wedge method).
     */
    std::string toString(CenterAlg center_alg);

    /**
     * \brief The way to aggregate the distances from curves to center within a
     * cluster (see paper).
     */
    enum class C2CDist {
        Median, Mean, Max
    };

    /**
     * \brief Compute the aggregated distance to the center curve for a cluster.
     * \param curves The set of curves clustered.
     * \param center_curve The center.
     * \param curve_ids The IDs of the curves in the cluster (subset of indices
     * of `curves').
     * \param c2c_dist The aggregation function.
     * \param dist The distance function for curves.
     * \return The result of computing the distances from each curve in the
     * cluster to the center, aggregated.
     */
    distance_t calcC2CDist(Curves const& curves, Curve const& center_curve,
        CurveIDs const& curve_ids, C2CDist c2c_dist,
        std::function<distance_t(Curve const&, Curve const&)> const& dist);
}

namespace detail {
    /**
     * \brief Update the centers. Best not to use directly.
     * \param curves The set of curves clustered.
     * \param clustering The clustering on curves (will be updated).
     * \param c2c_dist The aggregation function for clusters.
     * \param fix_start Whether to fix the starting points when updating the
     * center curves.
     * \param fix_end Ditto for the end points.
     * \param dist The distance function for the curves.
     * \param center_curve_update The single center curve update function.
     * \param args Extra arguments passed to the center update function.
     * \return True if centers have been updated.
     */
    template<typename F, typename... Args>
    bool updateCenters(clustering::Curves const& curves, Clustering& clustering,
            clustering::C2CDist c2c_dist, bool fix_start, bool fix_end,
            std::function<distance_t(Curve const&, Curve const&)> const& dist,
            F const& center_curve_update, Args&&... args) {

        bool found_new_center = false;
        for (auto& cluster: clustering)
            if (cluster.cost == std::numeric_limits<distance_t>::max())
                cluster.cost = calcC2CDist(curves, cluster.center_curve,
                    cluster.curve_ids, c2c_dist, dist);

        for (auto& cluster: clustering) {
            Curve new_center_curve = center_curve_update(curves, cluster,
                fix_start, fix_end, std::forward<Args>(args)...);

            if (cluster.center_curve != new_center_curve) {
                auto new_dist = calcC2CDist(curves, new_center_curve,
                    cluster.curve_ids, c2c_dist, dist);
                if (new_dist < cluster.cost) {
                    cluster.center_curve = std::move(new_center_curve);
                    cluster.cost = new_dist;
                    found_new_center = true;
                }
            }
        }
        return found_new_center;
    }
}

namespace clustering {
    /**
     * \brief Update the centers.
     * \param curves The set of curves clustered.
     * \param clustering The clustering on curves (will be updated).
     * \param center_alg The algorithm to use for center update.
     * \param fix_start Whether to fix the starting points when updating the
     * center curves.
     * \param fix_end Ditto for the end points.
     * \param dist The distance function for the curves.
     * \param args Extra arguments passed to the center update function.
     * \return True if centers have been updated.
     */
    template<typename... Args>
    bool computeCenters(Curves const& curves, Clustering& clustering,
            CenterAlg center_alg, bool fix_start, bool fix_end,
            std::function<distance_t(Curve const&, Curve const&)> const& dist,
            Args&&... args) {
        switch (center_alg) {
        case CenterAlg::fsa:
            return detail::updateCenters(curves, clustering, C2CDist::Max,
                fix_start, fix_end, dist, fsa_update);
        case CenterAlg::dba:
            return detail::updateCenters(curves, clustering, C2CDist::Mean,
                fix_start, fix_end, dist, dba_update);
        case CenterAlg::cdba:
            return detail::updateCenters(curves, clustering, C2CDist::Median,
                fix_start, fix_end, dist, cdba_update);
        case CenterAlg::wedge:
            return detail::updateCenters(curves, clustering, C2CDist::Median,
                fix_start, fix_end, dist, wedge_wrap<Args...>,
                std::forward<Args>(args)...);
        case CenterAlg::none:
            return false;
        default:
            ERROR("No matching center_alg enum passed.");
        }
    }
}
#endif
