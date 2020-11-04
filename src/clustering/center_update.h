#ifndef CENTER_UPDATE_H
#define CENTER_UPDATE_H

#include "basic_types.h"

namespace clustering {
    using Curves = std::vector<Curve>;

    /**
     * \brief Given a cluster and a set of curves, compute a new center for the
     * cluster with FSA (see the paper for details).
     * \param curves The set of curves.
     * \param cluster One of the clusters resulting from clustering curves.
     * \param fix_start Whether we should fix or are allowed to move the
     * starting point of the center curve.
     * \param fix_end Ditto for the end point.
     * \return The new center curve.
     */
    Curve fsa_update(Curves const& curves, Cluster const& cluster,
        bool fix_start = false, bool fix_end = false);

    /**
     * \brief Given a cluster and a set of curves, compute a new center for the
     * cluster with DBA (see the paper for details).
     * \param curves The set of curves.
     * \param cluster One of the clusters resulting from clustering curves.
     * \param fix_start Whether we should fix or are allowed to move the
     * starting point of the center curve.
     * \param fix_end Ditto for the end point.
     * \return The new center curve.
     */
    Curve dba_update(Curves const& curves, Cluster const& cluster,
        bool fix_start = false, bool fix_end = false);

    /**
     * \brief Given a cluster and a set of curves, compute a new center for the
     * cluster with CDBA (see the paper for details).
     * \param curves The set of curves.
     * \param cluster One of the clusters resulting from clustering curves.
     * \param fix_start Whether we should fix or are allowed to move the
     * starting point of the center curve.
     * \param fix_end Ditto for the end point.
     * \return The new center curve.
     */
    Curve cdba_update(Curves const& curves, Cluster const& cluster,
        bool fix_start = false, bool fix_end = false);

    /**
     * \brief Given a cluster and a set of curves, compute a new center for the
     * cluster with the Wedge method (see the paper for details).
     * \param curves The set of curves.
     * \param cluster One of the clusters resulting from clustering curves.
     * \param fix_start Whether we should fix or are allowed to move the
     * starting point of the center curve.
     * \param fix_end Ditto for the end point.
     * \param eps The step for the grid search for the midpoint of the wedge.
     * \param radius The number of steps taken in each direction for the grid
     * search.
     * \return The new center curve.
     */
    Curve wedge_update(Curves const& curves, Cluster const& cluster,
        bool fix_start = false, bool fix_end = false,
        distance_t eps = 0.125, int radius = 20);

    /**
     * \brief A wrapper to dispatch calls to wedge_update.
     *
     * If we try to use wedge_update directly, C++ does not know that it has
     * default arguments, so not setting the values for eps or radius leads
     * to code not compiling. This wrapper calls the function directly, so
     * default arguments work as expected.
     * \param curves The set of curves.
     * \param cluster One of the clusters resulting from clustering curves.
     * \param fix_start Whether we should fix or are allowed to move the
     * starting point of the center curve.
     * \param fix_end Ditto for the end point.
     * \param args The optional arguments for the wedge method: nothing, or just
     * eps, or eps and radius.
     * \return The new center curve.
     */
    template<typename... Args>
    Curve wedge_wrap(Curves const& curves, Cluster const& cluster,
            bool fix_start, bool fix_end, Args... args) {
        return wedge_update(curves, cluster, fix_start, fix_end,
            std::forward<Args>(args)...);
    }
}
#endif
