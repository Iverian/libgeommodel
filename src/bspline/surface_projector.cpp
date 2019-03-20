#include <bspline/distance_curve.hpp>
#include <bspline/distance_surface.hpp>
#include <bspline/surface_projector.hpp>
#include <gm/point.hpp>
#include <gm\surf_point.hpp>
#include <util/debug.hpp>

#include <limits>
#include <optional>
#include <stdexcept>

namespace gm {

SurfaceProjector::SurfaceProjector(const BSplineSurface::Impl& impl)
    : impl_(&impl)
    , patches_(impl.bezier_patches())
{
}

SurfPoint SurfaceProjector::call(const Point& p) const
{
    std::optional<SurfPoint> u = std::nullopt;
    auto d = std::numeric_limits<double>::max();

    for (auto i = std::begin(patches_);
         !cmp::zero(d) && i != std::end(patches_); ++i) {
        auto c = DistanceSurface(*i, p);

        while (true) {
            auto [umin, dmin] = c.min_init();
            d = std::min(d, dmin);

            if (c.is_candidate(d)) {
                auto v = bminimize(c, p, umin);
                if (v) {
                    auto fv
                        = sqr(dist(p, impl_->f(v.value()))); // c.f(v.value());
                    if (cmp::le(fv, d)) {
                        u = v;
                        d = fv;
                        break;
                    }
                    // break; ?
                }
                if (!c.eliminate_segment(d)) {
                    break;
                }
            } else {
                break;
            }
        }
    }

    check_if(u.has_value(), "Unable to project point {} on bspline surface",
             p);
    return u.value();
}

std::optional<SurfPoint> SurfaceProjector::bminimize(const DistanceSurface& c,
                                                     const Point& p,
                                                     const SurfPoint& r) const
    noexcept
{
    std::optional<SurfPoint> v;
    auto t = c.min_bord_check();

    if (t.first == WhereMin::FRONT && t.second == WhereMin::NIL) {
        v = minimize(p, {c.pfront().u, r.v}, c.pfront(), c.pback(), Mode::V);
    } else if (t.first == WhereMin::BACK && t.second == WhereMin::NIL) {
        v = minimize(p, {c.pback().u, r.v}, c.pfront(), c.pback(), Mode::V);
    } else if (t.first == WhereMin::NIL && t.second == WhereMin::FRONT) {
        v = minimize(p, {r.u, c.pfront().v}, c.pfront(), c.pback(), Mode::U);
    } else if (t.first == WhereMin::NIL && t.second == WhereMin::BACK) {
        v = minimize(p, {r.u, c.pback().v}, c.pfront(), c.pback(), Mode::U);
    } else if (t.first == WhereMin::FRONT && t.second == WhereMin::FRONT) {
        v = minimize(p, c.pfront(), c.pfront(), c.pback(), Mode::BOTH);
    } else if (t.first == WhereMin::FRONT && t.second == WhereMin::BACK) {
        v = minimize(p, {c.pfront().u, c.pback().v}, c.pfront(), c.pback(),
                     Mode::BOTH);
    } else if (t.first == WhereMin::BACK && t.second == WhereMin::FRONT) {
        v = minimize(p, {c.pback().u, c.pfront().v}, c.pfront(), c.pback(),
                     Mode::BOTH);
    } else if (t.first == WhereMin::BACK && t.second == WhereMin::BACK) {
        v = minimize(p, c.pback(), c.pfront(), c.pback(), Mode::BOTH);
    } else {
        v = minimize(p, r, c.pfront(), c.pback(), Mode::BOTH);
    }

    return v;
}

std::optional<SurfPoint>
SurfaceProjector::minimize(const Point& p, SurfPoint r, const SurfPoint& a,
                           const SurfPoint& b, Mode m) const noexcept
{
    static constexpr auto max_iter = size_t(20);

    for (size_t i = 0; i < max_iter; ++i) {
        auto w = Vec(p, impl_->f(r));
        auto u = impl_->dfu(r);
        auto v = impl_->dfv(r);
        if (cmp::zero(w) || cmp::zero(sin(w, cross(u, v)))) {
            return r;
        }

        auto h = bord_check(r + next_step(r, w, u, v, m), a, b) - r;
        if (cmp::zero(u * h.u + v * h.v)) {
            break;
        }
        r += h;
    }

    return std::nullopt;
}

SurfPoint SurfaceProjector::next_step(SurfPoint r, const Vec& w, const Vec& fu,
                                      const Vec& fv, Mode m) const noexcept
{
    SurfPoint result;

    switch (m) {
    case Mode::BOTH: {
        double a[]
            = {dot(impl_->dfuu(r), w) + sqr(fu),
               dot(impl_->dfuv(r), w) + dot(fu, fv),
               dot(impl_->dfvv(r), w) + sqr(fv), -dot(fu, w), -dot(fv, w)};
        double d[] = {a[0] * a[2] - sqr(a[1]), a[3] * a[2] - a[4] * a[1],
                      a[4] * a[0] - a[3] * a[1]};
        result = {d[1] / d[0], d[2] / d[0]};
        break;
    }
    case Mode::U: {
        result = {-dot(fu, w) / (dot(impl_->dfuu(r), w) + sqr(fu)), 0};
        break;
    }
    case Mode::V: {
        result = {0, -dot(fv, w) / (dot(impl_->dfvv(r), w) + sqr(fv))};
        break;
    }
    }

    return result;
}

SurfPoint SurfaceProjector::bord_check(SurfPoint r, const SurfPoint& a,
                                       const SurfPoint& b) const noexcept
{
    if (r.u < a.u) {
        r.u = a.u;
    } else if (r.u > b.u) {
        r.u = b.u;
    }
    if (r.v < a.v) {
        r.v = a.v;
    } else if (r.v > b.v) {
        r.v = b.v;
    }
    return r;
}

} // namespace gm