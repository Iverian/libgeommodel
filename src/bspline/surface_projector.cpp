#define _USE_MATH_DEFINES

#include <bspline/distance_curve.hpp>
#include <bspline/distance_surface.hpp>
#include <bspline/surface_projector.hpp>
#include <bspline/util.hpp>
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
    auto d = init_global(p);

    for (auto i = std::begin(patches_);
         !cmp::zero(d) && i != std::end(patches_); ++i) {
        auto c = DistanceSurface(*i, p);

        auto [uloc, dloc] = min_init(p, c);
        if (dloc < d) {
            u = uloc;
            d = dloc;
        }
        while (c.is_candidate(d)) {
            auto [init, mode] = init_value(c, uloc);
            auto v = minimize(c, p, init, mode);
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

std::pair<SurfPoint, SurfaceProjector::OptMode>
SurfaceProjector::init_value(const DistanceSurface& c, SurfPoint r) const
    noexcept
{
    OptMode m = OptMode::BOTH;
    auto b = c.is_min_on_border();

    if (b) {
        switch (b.value()) {
        case WhereMin::UFRONT: {
            r.u = c.pfront().u;
            m = OptMode::V;
            break;
        }
        case WhereMin::UBACK: {
            r.u = c.pback().u;
            m = OptMode::V;
            break;
        }
        case WhereMin::VFRONT: {
            r.v = c.pfront().v;
            m = OptMode::U;
            break;
        }
        case WhereMin::VBACK: {
            r.v = c.pback().v;
            m = OptMode::U;
            break;
        }
        }
    }

    return {r, m};
}

std::optional<SurfPoint> SurfaceProjector::minimize(const DistanceSurface& c,
                                                    const Point& p,
                                                    SurfPoint r,
                                                    OptMode m) const noexcept
{
    static constexpr auto max_iter = size_t(25);

    std::optional<SurfPoint> result = std::nullopt;
    auto a = c.pfront();
    auto b = c.pback();

    Vec w, u, v;
    SurfPoint s;
    for (size_t i = 0; i < max_iter; ++i) {
        w = Vec(p, impl_->f(r));
        u = impl_->dfu(r);
        v = impl_->dfv(r);
        if (cmp::zero(w) || (cmp::zero(cos(u, w)) && cmp::zero(cos(v, w)))) {
            result = r;
            break;
        }

        s = bord_check(r + next_step(r, w, u, v, m), a, b);
        if (cmp::zero((s.u - r.u) * u + (s.v - r.v) * v)) {
            result = r;
            break;
        }
        r = s;
    }

    return result;
}

SurfPoint SurfaceProjector::next_step(SurfPoint r, const Vec& w, const Vec& fu,
                                      const Vec& fv, OptMode m) const noexcept
{
    SurfPoint result;

    switch (m) {
    case OptMode::BOTH: {
        double a[]
            = {dot(impl_->dfuu(r), w) + sqr(fu),
               dot(impl_->dfuv(r), w) + dot(fu, fv),
               dot(impl_->dfvv(r), w) + sqr(fv), -dot(fu, w), -dot(fv, w)};
        double d[] = {a[0] * a[2] - sqr(a[1]), a[3] * a[2] - a[4] * a[1],
                      a[4] * a[0] - a[3] * a[1]};
        result = {d[1] / d[0], d[2] / d[0]};
        break;
    }
    case OptMode::U: {
        result = {-dot(fu, w) / (dot(impl_->dfuu(r), w) + sqr(fu)), 0};
        break;
    }
    case OptMode::V: {
        result = {0, -dot(fv, w) / (dot(impl_->dfvv(r), w) + sqr(fv))};
        break;
    }
    }

    return result;
}

SurfPoint SurfaceProjector::bord_check(SurfPoint r, const SurfPoint& a,
                                       const SurfPoint& b) const noexcept
{
    r.u = ::bord_check(r.u, a.u, b.u);
    r.v = ::bord_check(r.v, a.v, b.v);
    return r;
}

double SurfaceProjector::rf(const Point& p, const SurfPoint& r) const noexcept
{
    return sqr(Vec(p, impl_->f(r)));
}

double SurfaceProjector::init_global(const Point& p) const
{
    auto result = std::numeric_limits<double>::max();
    auto& ku = impl_->knots().first;
    auto& kv = impl_->knots().second;

    for (auto u : {ku.front(), ku.back()}) {
        for (auto v : {kv.front(), kv.back()}) {
            auto d = rf(p, {u, v});
            if (result < d) {
                result = d;
            }
        }
    }

    return result;
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

#define norm(x) (cross(impl_->dfu(x), impl_->dfv(x)))

bool SurfaceProjector::is_flat(const DistanceSurface& c) const noexcept
{
    static constexpr auto upper_angle = M_PI / 8;

    auto s = c.shape();
    auto m = norm((c.pfront() + c.pback()) / 2);
    auto a = std::numeric_limits<double>::max();

    for (size_t i = 0; i < s.first; ++i) {
        for (size_t j = 0; j < s.second; ++j) {
            auto r = c.itarg({i, j});
            auto n = norm(r);
            if (auto cur = cos(m, n); cur < a) {
                a = cur;
            }
        }
    }

    return cmp::le(acos(a), upper_angle);
}

#undef norm

} // namespace gm