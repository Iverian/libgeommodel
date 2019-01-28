#include "curve_projector.h"

#include <util/debug.h>
#include <util/math.h>

#include <algorithm>

using namespace std;

namespace gm {

double binom(size_t n, size_t k) __GM_NOEXCEPT_RELEASE__;

template <class T>
const T& get(const vector<T>& vec, size_t index)
{
    static const T def {};

    return (index < vec.size()) ? vec[index] : def;
}

CurveProjector::CurveProjector(const BSplineCurve::Impl& impl)
    : patches_(impl.bezier_patches())
    , front_(patches_.front().cp.front().p())
    , back_(patches_.back().cp.back().p())
{
}

double CurveProjector::get(const Point& p) const
{
    ::vector<DistanceCurve> dcurves(patches_.size());
    ::transform(::begin(patches_), ::end(patches_), ::begin(dcurves),
                [&p](auto& c) { return DistanceCurve(c, p); });

    double mdist, result;
    {
        auto fd = dist(p, front_);
        auto bd = dist(p, back_);
        if (fd < bd || isnear(fd, bd)) {
            mdist = fd;
            result = patches_.front().front;
        } else {
            mdist = bd;
            result = patches_.back().back;
        }
    }

    for (auto& c : dcurves) {
        if (is_candidate(mdist, c)) {
            if (c.peak_point().has_value()) {
                auto it = c.minimize((c.pfront() + c.pback()) / 2);
                if (it.has_value()) {
                    result = it.value();
                    mdist = c.f(result);
                }
            } else {
            }
        }
    }

    throw runtime_error("not implemented");
    return 0;
}

bool CurveProjector::is_candidate(double a, const DistanceCurve& c) const
    noexcept
{
    return ::any_of(::begin(c.points()), ::end(c.points()), [&a](auto& wp) {
        auto p = wp.p()[0];
        return p < a || isnear(p, a);
    });
}

DistanceCurve::DistanceCurve()
    : order_(0)
    , knots_()
    , points_()
    , cdb_()
{
}

DistanceCurve::DistanceCurve(const BezierPatch& patch, const Point& point)
    : order_(2 * (patch.order() - 1) + 1)
    , knots_(2 * order_)
    , points_(order_)
    , cdb_()
{
    auto pr = point.raw();
    auto deg = patch.order() - 1;
    double c = 0;

    for (size_t i = 0; i < order_; ++i) {
        knots_[i] = patch.front;
        knots_[order_ + i] = patch.back;

        auto& q = points_[i];
        auto first = (i < deg) ? size_t(0) : (i - deg);
        auto last = ::min(deg, i) + 1;
        for (auto k = first; k < last; ++k) {
            auto a = patch.cp[k] - pr;
            auto b = patch.cp[i - k] - pr;

            c = binom(deg, k) * binom(deg, i - k);
            q += c * wdot(a, b);
        }
        c = binom(2 * deg, i);
        q /= c;
    }

    cdb_ = CoxDeBoor<DPoint>(order_, knots_, points_);
}

DistanceCurve& DistanceCurve::reconfig(DistanceCurve& out,
                                       double a) const __GM_NOEXCEPT_RELEASE__
{
    check_ifd(out.order_ == order_, "Output curve is incompatible");
    ::copy(::begin(knots_), ::end(knots_), begin(out.knots_));
    ::copy(::begin(points_), ::end(points_), ::begin(out.points_));

    for (auto& p : out.points_) {
        p[0] -= p[1] * a;
    }
    return out;
}

size_t DistanceCurve::order() const noexcept
{
    return order_;
}

::vector<double> DistanceCurve::knots() const noexcept
{
    return knots_;
}

::vector<DPoint> DistanceCurve::points() const noexcept
{
    return points_;
}

double DistanceCurve::pfront() const noexcept
{
    return knots_.front();
}

double DistanceCurve::pback() const noexcept
{
    return knots_.back();
}

double DistanceCurve::f(double u) const noexcept
{
    return cdb_.proxy(u).get(0).p()[0];
}

double DistanceCurve::df(double u) const noexcept
{
    auto p = cdb_.proxy(u).range(2);
    return DPoint::d1(p).p()[0];
}

double DistanceCurve::df2(double u) const noexcept
{
    auto p = cdb_.proxy(u).range(3);
    return DPoint::d2(p).p()[0];
}

bool DistanceCurve::min_between(double a, double b) const noexcept
{
    return ::signbit(df(a)) != ::signbit(df(b));
}

::optional<double> DistanceCurve::minimize(double u0) const noexcept
{
    static constexpr size_t max_iter = 10000;

    auto u = u0;
    for (size_t i = 0; i < max_iter; ++i) {
        auto h = bound_check(u - df(u) / df2(u)) - u;
        auto v = armijo_step(u, h);

        if (isnear(u, v, Tolerance::MAX)) {
            return v;
        }
        u = v;
    }

    return nullopt;
}

::optional<size_t> DistanceCurve::peak_point() const noexcept
{
    bool flag = false;
    size_t result = 0;

    double prev, cur = points_[0].p()[0];
    for (size_t i = 1; i < order_; ++i) {
        prev = cur;
        cur = points_[i].p()[0];

        if (prev < cur) {
            if (flag) {
                return nullopt;
            }
        } else if (!isnear(prev, cur)) {
            if (!flag) {
                result = i;
                flag = true;
            }
        } else {
            return nullopt;
        }
    }

    return result;
}

double DistanceCurve::bound_check(double u) const noexcept
{
    if (u > knots_.back()) {
        u = knots_.back();
    } else if (u < knots_.front()) {
        u = knots_.front();
    }

    return u;
}

double DistanceCurve::armijo_step(double u, double h) const noexcept
{
    static constexpr size_t max_iter = 10000;
    static constexpr double a = 0.5;
    static constexpr double b = 0.9;

    auto fu = f(u);
    auto d = a * h * df(u);
    double k = 1.;

    for (size_t i = 0; i < max_iter; ++i) {
        if (f(u + b * k * h) - fu > b * k * d) {
            return u + k * h;
        }
        k *= b;
    }

    return u + h;
}

inline double binom(size_t n, size_t k) __GM_NOEXCEPT_RELEASE__
{
    check_ifd(k <= n, "Unable to compute C^n_k with n = {}, k = {}", n, k);

    double result = 1;
    if (k != 0 && k != n) {
        for (size_t i = 0; i < k; ++i) {
            result *= double(n - i) / (k - i);
        }
    }
    return result;
}

} // namespace gm