#define _USE_MATH_DEFINES

#include <bspline/distance_curve.hpp>
#include <bspline/distance_surface.hpp>
#include <bspline/surface_projector.hpp>
#include <gm/compare.hpp>
#include <gm/point.hpp>
#include <gm/surf_point.hpp>
#include <gm/vec.hpp>
#include <util/cyclic_iterator.hpp>
#include <util/debug.hpp>

#include <algorithm>
#include <cmath>
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

        auto [uloc, dloc] = min_init(p, c);
        if (dloc < d) {
            u = uloc;
            d = dloc;
        }
        while (c.is_candidate(d)) {
            auto v = minimize(c, p, uloc);
            if (v) {
                if (auto fv = rf(p, v.value()); cmp::le(fv, d)) {
                    u = v;
                    d = fv;
                    break;
                }
            }
            if (!c.eliminate_segment(d)) {
                break;
            }

            std::tie(uloc, dloc) = min_init(p, c);
            if (dloc < d) {
                u = uloc;
                d = dloc;
            }
        }
    }

    check_if(u, "Unable to project point {} on bspline surface", p);
    return u.value();
}

std::optional<SurfPoint> SurfaceProjector::minimize(const DistanceSurface& c,
                                                    const Point& p,
                                                    SurfPoint r) const noexcept
{
    static constexpr auto max_iter = size_t(100);

    std::optional<SurfPoint> result;
    Vec w, u, v;
    SurfPoint s;
    auto a = c.pfront();
    auto b = c.pback();

    for (size_t i = 0; i < max_iter; ++i) {
        w = Vec(p, impl_->f(r));
        u = impl_->dfu(r);
        v = impl_->dfv(r);
        if (cmp::zero(w) || cmp::zero(angle(w, cross(u, v)))) {
            result = r;
            break;
        }

        s = bord_check(r + next_step(r, w, u, v), a, b);
        if (cmp::zero(w) || (cmp::zero(cos(u, w)) && cmp::zero(cos(v, w)))
            || cmp::zero((s.u - r.u) * u + (s.v - r.v) * v)) {
            result = r;
            break;
        }
        r = s;
    }

    return result;
}

SurfPoint SurfaceProjector::next_step(SurfPoint r, const Vec& w, const Vec& fu,
                                      const Vec& fv) const noexcept
{
    double a[] = {dot(impl_->dfuu(r), w) + sqr(fu),
                  dot(impl_->dfuv(r), w) + dot(fu, fv),
                  dot(impl_->dfvv(r), w) + sqr(fv), -dot(fu, w), -dot(fv, w)};
    double d[] = {a[0] * a[2] - sqr(a[1]), a[3] * a[2] - a[4] * a[1],
                  a[4] * a[0] - a[3] * a[1]};
    return {d[1] / d[0], d[2] / d[0]};
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

double SurfaceProjector::rf(const Point& p, const SurfPoint& r) const noexcept
{
    return sqr(impl_->f(r) - p);
}

std::pair<SurfPoint, double>
SurfaceProjector::min_init(const Point& p, const DistanceSurface& c) const
    noexcept
{
    auto s = c.shape();
    auto dmin = std::numeric_limits<double>::max();
    SurfPoint argmin {};

    for (size_t i = 0; i < s.first; ++i) {
        for (size_t j = 0; j < s.second; ++j) {
            auto r = c.itarg({i, j});
            auto d = rf(p, r);
            if (d < dmin) {
                argmin = r;
                dmin = d;
            }
        }
    }

    return {argmin, dmin};
}

} // namespace gm