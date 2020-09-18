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
     * \param fix_endpoints Whether we should fix or are allowed to move the
     * endpoints of the center curve.
     * \return The new center curve.
     */
    Curve fsa_update(Curves const& curves, Cluster const& cluster,
        bool fix_endpoints = false);

    /**
     * \brief Given a cluster and a set of curves, compute a new center for the
     * cluster with DBA (see the paper for details).
     * \param curves The set of curves.
     * \param cluster One of the clusters resulting from clustering curves.
     * \param fix_endpoints Whether we should fix or are allowed to move the
     * endpoints of the center curve.
     * \return The new center curve.
     */
    Curve dba_update(Curves const& curves, Cluster const& cluster,
        bool fix_endpoints = false);

    /**
     * \brief Given a cluster and a set of curves, compute a new center for the
     * cluster with CDBA (see the paper for details).
     * \param curves The set of curves.
     * \param cluster One of the clusters resulting from clustering curves.
     * \param fix_endpoints Whether we should fix or are allowed to move the
     * endpoints of the center curve.
     * \return The new center curve.
     */
    Curve cdba_update(Curves const& curves, Cluster const& cluster,
        bool fix_endpoints = false);

    /**
     * \brief Given a cluster and a set of curves, compute a new center for the
     * cluster with the Wedge method (see the paper for details).
     * \param curves The set of curves.
     * \param cluster One of the clusters resulting from clustering curves.
     * \param fix_endpoints Whether we should fix or are allowed to move the
     * endpoints of the center curve.
     * \return The new center curve.
     */
    Curve wedge_update(Curves const& curves, Cluster const& cluster,
        bool fix_endpoints = false, distance_t eps = 0.125, int radius = 20);
}
#endif
