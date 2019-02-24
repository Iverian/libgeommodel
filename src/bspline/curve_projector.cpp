#include <bspline/array_operators.h>
#include <bspline/curve_projector.h>
#include <bspline/wpoint.h>
#include <gm/dot.h>
#include <util/cyclic_iterator.h>
#include <util/debug.h>
#include <util/math.h>

#include <algorithm>
#include <array>
#include <iterator>
#include <vector>

using namespace std;

namespace gm {

CurveProjector::CurveProjector(const BSplineCurve::Impl::Super& impl)
    : parent_(&impl)
    , patches_(impl.bezier_patches())
{
}

double CurveProjector::call(const Point& p) const
{
    ::vector<DistanceCurve> dcurves(patches_.size());
    ::vector<bool> eliminated(dcurves.size());
    ::optional<double> u = nullopt;
    auto d = ::numeric_limits<double>::max();
    auto r = p.raw();

    ::transform(::begin(patches_), ::end(patches_), ::begin(dcurves),
                [&p](auto& c) { return DistanceCurve(c, p); });

    auto it = CyclicIterator(::begin(dcurves), ::end(dcurves));
    auto to_eliminate = dcurves.size();

    for (; to_eliminate != 0; ++it) {
        auto b = eliminated[::distance(it.first(), it.iter())];
        if (!b) {
            auto [umin, dmin] = it->min_init();
            d = ::min(d, dmin);

            if (it->is_candidate(d)) {
                if (it->peak_point()) {
                    if (auto v = minimize(r, umin, it->pfront(), it->pback());
                        v) {
                        if (auto fv = it->f(v.value()); cmp::le(fv, d)) {
                            u = v;
                            d = fv;
                            if (cmp::zero(d)) {
                                break;
                            }
                        }
                        --to_eliminate;
                        b = true;

                        continue;
                    }
                }
                if (!it->eliminate_segment(d)) {
                    --to_eliminate;
                    b = true;
                }
            } else {
                --to_eliminate;
                b = true;
            }
        }
    }

    return u.value();
}

std::optional<double>
CurveProjector::minimize(const BSplineCurve::Impl::Super::CPoint::Proj& p,
                         double u0, const double& pfront,
                         const double& pback) const noexcept
{
    static constexpr auto max_iter = size_t(25);

    auto u = u0;
    for (size_t i = 0; i < max_iter; ++i) {
        auto w = Vec(parent_->f(u) - p);
        auto d = Vec(parent_->df(u));
        auto d2 = Vec(parent_->df2(u));
        if (cmp::zero(w) || cmp::zero(cos(d, w))) {
            return u;
        }

        auto h
            = bord_check(u - dot(d, w) / (dot(d2, w) + sqr(d)), pfront, pback)
            - u;
        if (cmp::zero(h * d)) {
            break;
        }

        u += h;
    }

    return nullopt;
}

double CurveProjector::bord_check(double u, const double& a,
                                  const double& b) const noexcept
{
    while (u < a || u > b) {
        if (u < a) {
            u = 2 * a - u;
        }
        if (u > b) {
            u = 2 * b - u;
        }
    }
    return u;
}

double CurveProjector::armijo_step(
    double u, double h, const Vec& w, const Vec& d,
    const BSplineCurve::Impl::Super::CPoint::Proj& p) const noexcept
{
    static constexpr auto a = 0.5;
    static constexpr auto b = 0.9;

    auto fu = sqr(w);
    auto rhs = 2 * a * h * dot(w, d);
    double j = 1.;
    while (true) {
        auto k = b * j;
        auto z = Vec(parent_->f(u + k * h) - p);
        if (sqr(z) > k * rhs + fu) {
            break;
        }
        j = k;
    }
    return j;
}

} // namespace gm