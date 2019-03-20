#include <bspline/curve_projector.hpp>
#include <bspline/wpoint.hpp>
#include <gm/dot.hpp>
#include <gm/point.hpp>
#include <util/debug.hpp>
#include <util/math.hpp>

#include <algorithm>
#include <array>
#include <iterator>
#include <vector>

namespace gm {

CurveProjector::CurveProjector(const BSplineCurve::Impl& impl)
    : impl_(&impl)
    , patches_(impl.bezier_patches())
{
}

double CurveProjector::call(const Point& p) const
{
    std::optional<double> u = std::nullopt;
    auto d = std::numeric_limits<double>::max();

    for (auto i = std::begin(patches_);
         i != std::end(patches_) && !cmp::zero(d); ++i) {
        auto c = DistanceCurve(*i, p);

        while (true) {
            auto [umin, dmin] = c.min_init();
            d = std::min(d, dmin);

            if (c.is_candidate(d)) {
                if (c.peak_point()) {
                    auto v = minimize(p, umin, c.pfront(), c.pback());
                    if (v.has_value()) {
                        auto fv = c.f(v.value());
                        if (cmp::le(fv, d)) {
                            u = v;
                            d = fv;
                        }
                        break;
                    }
                }
                if (!c.eliminate_segment(d)) {
                    break;
                }
            } else {
                break;
            }
        }
    }

    check_if(u.has_value(), "Unable to project point {} on bspline curve", p);
    return u.value();
}

std::optional<double> CurveProjector::minimize(const Point& p, double u0,
                                               const double& a,
                                               const double& b) const noexcept
{
    static constexpr auto max_iter = size_t(20);

    double u = u0, h = 0.;
    for (size_t i = 0; i < max_iter; ++i) {
        auto w = Vec(p, impl_->f(u));
        auto d = impl_->df(u);
        auto g = impl_->df2(u);
        if (cmp::zero(w) || cmp::zero(cos(d, w))) {
            return u;
        }

        h = bord_check(u - dot(d, w) / (dot(g, w) + sqr(d)), a, b) - u;
        if (cmp::zero(h * d)) {
            break;
        }

        u += h;
    }

    return std::nullopt;
}

double CurveProjector::bord_check(double u, const double& a,
                                  const double& b) const noexcept
{
    if (u < a) {
        u = a;
    }
    if (u > b) {
        u = b;
    }
    return u;
}

} // namespace gm