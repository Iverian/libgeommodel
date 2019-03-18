#include <bspline/distance_curve.hpp>
#include <bspline/distance_surface.hpp>
#include <bspline/surface_projector.hpp>
#include <gm/point.hpp>

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
                auto v = minimize(p, umin, c.pfront(), c.pback());
                if (v) {
                    auto fv = c.f(v.value());
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

    throw std::runtime_error("not implemented");
    return SurfPoint();
}

std::optional<SurfPoint>
SurfaceProjector::minimize(const Point& p, SurfPoint r, const SurfPoint& a,
                           const SurfPoint& b) const noexcept
{
    static constexpr auto max_iter = size_t(100);

    for (size_t i = 0; i < max_iter; ++i) {
        auto w = Vec(p, impl_->f(r));
        auto u = impl_->dfu(r);
        auto v = impl_->dfv(r);
        if (cmp::zero(w) || cmp::zero(sin(w, cross(u, v)))) {
            return r;
        }

        auto h = bord_check(r + next_step(r, w, u, v), a, b) - r;
        if (cmp::zero(u * h.u + v * h.v)) {
            break;
        }
        r += h;
    }

    return std::nullopt;
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

} // namespace gm