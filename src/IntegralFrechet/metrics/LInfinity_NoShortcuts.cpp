#include "IntegralFrechet/metrics/include.h"

namespace {
    void join_paths(Points& path_forward, const Points& path_backward) {
        if (!approx_equal(path_forward.back(), path_backward.back())) {
            path_forward.push_back(path_backward.back());
        }
        path_forward.insert(path_forward.end(), path_backward.rbegin() + 1, path_backward.rend());
    }

    void join_forward_paths(Points& path_left, const Points& path_right) {
        assert(approx_equal(path_left.back(), path_right.front()));
        path_left.insert(path_left.end(), path_right.begin() + 1, path_right.end());
    }

    Point linf_sd_diag(Points& path, const Cell& cell, const Point& s, const Point& t, BFDirection dir) {
        const auto line_left = (dir == BFDirection::Forward ? -1 : 1);
        if (cell.ell_v.side(s) != -line_left || cell.ell_h.side(s) != line_left) {
            return s;
        }

        MonotoneComparator compare(dir);
        const Line diag(s, {1, 1});
        auto new_s = std::min({
            intersect(diag, Line::vertical(t)),
            intersect(diag, Line::horizontal(t)),
            intersect(diag, cell.ell_h),
            intersect(diag, cell.ell_v),
        }, compare);

        if (!approx_equal(new_s, s)) {
            path.push_back(new_s);
            return new_s;
        }
        return s;
    }

    Point linf_sd_hv(Points& path, const Cell& cell, const Point& s, const Point& t, BFDirection dir) {
        const auto line_left = (dir == BFDirection::Forward ? -1 : 1);
        MonotoneComparator compare(dir);

        Point new_s = t;

        if (cell.ell_v.side(s) != -line_left && cell.ell_h.side(s) == line_left) {
            if (compare(s, cell.mid)) {
                if (!cell.ell_v.isHorizontal()) {
                    new_s = std::min(new_s, intersect(Line::horizontal(s), cell.ell_v), compare);
                }
            } else {
                new_s = std::min(new_s, intersect(Line::horizontal(s), cell.ell_h), compare);
            }
            new_s = std::min(new_s, intersect(Line::horizontal(s), Line::vertical(t)), compare);
        } else if (cell.ell_h.side(s) != line_left && cell.ell_v.side(s) == -line_left) {
            if (compare(s, cell.mid)) {
                if (!cell.ell_h.isVertical()) {
                    new_s = std::min(new_s, intersect(Line::vertical(s), cell.ell_h), compare);
                }
            } else {
                new_s = std::min(new_s, intersect(Line::vertical(s), cell.ell_v), compare);
            }
            new_s = std::min(new_s, intersect(Line::vertical(s), Line::horizontal(t)), compare);
        } else {
            return s;
        }

        if (!approx_equal(new_s, s)) {
            path.push_back(new_s);
            return new_s;
        }
        return s;
    }

    Point linf_sd_ell(Points& path, const Cell& cell, const Point& s, const Point& t, BFDirection dir) {
        MonotoneComparator compare(dir);

        if (!compare(s, cell.mid)) {
            return s;
        }

        Point new_s;

        if (cell.ell_v.side(s) == 0) {
            new_s = std::min({
                cell.mid,
                intersect(cell.ell_v, Line::vertical(t)),
            }, compare);

            if (!cell.ell_v.isHorizontal()) {
                new_s = std::min(new_s, intersect(cell.ell_v, Line::horizontal(t)), compare);
            }
        } else if (cell.ell_h.side(s) == 0) {
            new_s = std::min({
                cell.mid,
                intersect(cell.ell_h, Line::horizontal(t)),
            }, compare);

            if (!cell.ell_h.isVertical()) {
                new_s = std::min(new_s, intersect(cell.ell_h, Line::vertical(t)), compare);
            }
        } else {
            return s;
        }

        if (!approx_equal(new_s, s)) {
            path.push_back(new_s);
            return new_s;
        }
        return s;
    }
}

template<>
Points compute_matching<ParamMetric::LInfinity_NoShortcuts>(const Cell& cell, const Point& cell_s, const Point& cell_t) {
    const auto& ell_h = cell.ell_h;
    const auto& ell_v = cell.ell_v;

    const bool is_opposite_direction = !approx_zero(cell.ell_v.direction.y) && cell.ell_v.direction.y < 0;

    const auto cell_len1 = std::abs(cell_s.x - cell_t.x);
    const auto cell_len2 = std::abs(cell_s.y - cell_t.y);

    auto s = cell_s;
    auto t = cell_t;

    Points path_forward = {s};
    Points path_backward = {t};

    // Return if trivial
    if (approx_equal(s.x, t.x) || approx_equal(s.y, t.y)) {
        join_paths(path_forward, path_backward);
        return path_forward;
    }

    assert(cell_s.x <= cell_t.x);
    assert(cell_s.y <= cell_t.y);
    MonotoneComparator compare(BFDirection::Forward);

    // Check if we need to deal with the case where the steepest descent
    // paths would not meet in a single point:
    if (is_opposite_direction && !approx_equal(cell_len1, cell_len2)) {
        if (cell_len2 > cell_len1) {
            if (ell_h.side(cell_s) != 0 && ell_h.side(cell_t) != 0 && ell_h.side(cell_s) != ell_h.side(cell_t)) {
                // Vertical plateau case

                const Line diag(
                    cell_s + Point(0, (cell_len2 - cell_len1) / 2),
                    {1, 1}
                );

                auto p = intersect(diag, ell_h);

                if (!compare(cell_s, p)) {
                    p = std::max({
                        intersect(diag, Line::horizontal(cell_s)),
                        intersect(diag, Line::vertical(cell_s)),
                    }, compare);
                } else if (!compare(p, cell_t)) {
                    p = std::max({
                        intersect(diag, Line::horizontal(cell_t)),
                        intersect(diag, Line::vertical(cell_t)),
                    }, compare);
                }

                Points path_left = compute_matching<ParamMetric::LInfinity_NoShortcuts>(cell, cell_s, p);
                Points path_right = compute_matching<ParamMetric::LInfinity_NoShortcuts>(cell, p, cell_t);

                join_forward_paths(path_left, path_right);
                return path_left;
            }
        } else {
            if (ell_v.side(cell_s) != 0 && ell_v.side(cell_t) != 0 && ell_v.side(cell_s) != ell_v.side(cell_t)) {
                // Horizontal plateau case

                const Line diag(
                    cell_s + Point((cell_len1 - cell_len2) / 2, 0),
                    {1, 1}
                );

                auto p = intersect(diag, ell_v);

                if (!compare(cell_s, p)) {
                    p = std::max({
                        intersect(diag, Line::horizontal(cell_s)),
                        intersect(diag, Line::vertical(cell_s)),
                    }, compare);
                } else if (!compare(p, cell_t)) {
                    p = std::max({
                        intersect(diag, Line::horizontal(cell_t)),
                        intersect(diag, Line::vertical(cell_t)),
                    }, compare);
                }

                Points path_left = compute_matching<ParamMetric::LInfinity_NoShortcuts>(cell, cell_s, p);
                Points path_right = compute_matching<ParamMetric::LInfinity_NoShortcuts>(cell, p, cell_t);

                join_forward_paths(path_left, path_right);
                return path_left;
            }
        }
    }

    // Try steepest descent parallel to middle ellipse axis
    s = linf_sd_diag(path_forward, cell, s, t, BFDirection::Forward);
    t = linf_sd_diag(path_backward, cell, t, s, BFDirection::Backward);

    // Return if trivial
    if (approx_equal(s.x, t.x) || approx_equal(s.y, t.y)) {
        join_paths(path_forward, path_backward);
        return path_forward;
    }

    // Try steepest descent parallel to cell axes
    s = linf_sd_hv(path_forward, cell, s, t, BFDirection::Forward);
    t = linf_sd_hv(path_backward, cell, t, s, BFDirection::Backward);

    // Return if trivial
    if (approx_equal(s.x, t.x) || approx_equal(s.y, t.y)) {
        join_paths(path_forward, path_backward);
        return path_forward;
    }

    // Try steepest descent along other ellipse axes
    s = linf_sd_ell(path_forward, cell, s, t, BFDirection::Forward);
    t = linf_sd_ell(path_backward, cell, t, s, BFDirection::Backward);

    // Return if trivial
    if (approx_equal(s.x, t.x) || approx_equal(s.y, t.y)) {
        join_paths(path_forward, path_backward);
        return path_forward;
    }

    throw std::logic_error("should have found a path by now");
}

template<>
distance_t integrate_linear_dist<ParamMetric::LInfinity_NoShortcuts>(const Cell& /*cell*/, const Point& s, const Point& t) {
    return std::max(std::abs(t.x - s.x), std::abs(t.y - s.y));
}
