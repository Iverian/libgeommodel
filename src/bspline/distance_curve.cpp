#include <bspline/distance_curve.hpp>
#include <bspline/util.hpp>
#include <bspline/wpoint.hpp>
#include <util/cyclic_iterator.hpp>
#include <util/math.hpp>
#include <gm/misc.hpp>

#include <iterator>
#include <optional>

namespace gm {

DistanceCurve::DistanceCurve()
    : c_()
{
}

DistanceCurve::DistanceCurve(const BSplineCurve::Impl::BezierPatch& patch,
                             const Point& point)
    : c_()
{
    auto order = 2 * patch.order() - 1;
    Super::KnotsType knots(2 * order);
    Super::CPointsType cp(order);

    auto pr = point.raw();
    auto deg = patch.order() - 1;
    double c = 0;

    for (size_t i = 0; i < order; ++i) {
        knots[i] = patch.pfront();
        knots[order + i] = patch.pback();

        auto& q = cp[i];
        auto first = (i < deg) ? size_t(0) : (i - deg);
        auto last = std::min(deg, i) + 1;
        for (auto k = first; k < last; ++k) {
            auto a = patch.cpoints()[k] - pr;
            auto b = patch.cpoints()[i - k] - pr;

            c = binom(deg, k) * binom(deg, i - k);
            q += c * wdot(a, b);
        }
        c = binom(2 * deg, i);
        q /= c;
    }

    c_ = Super(order, std::move(knots), std::move(cp));
}

double DistanceCurve::pfront() const noexcept
{
    return c_.pfront();
}

double DistanceCurve::pback() const noexcept
{
    return c_.pback();
}

double DistanceCurve::itarg(size_t cp_index) const noexcept
{
    return pfront() + cp_index * (pback() - pfront()) / (c_.order() - 1);
}

double DistanceCurve::argti(double arg) const noexcept
{
    return (c_.order() - 1) * (arg - pfront()) / (pback() - pfront());
}

double DistanceCurve::tocparg(double arg, bool dir) const noexcept
{
    arg = argti(arg);
    return itarg(dir ? ceil<size_t>(arg) : floor<size_t>(arg));
}

double DistanceCurve::f(double u) const noexcept
{
    return c_.f(u)[0];
}

std::pair<double, double> DistanceCurve::min_init() const noexcept
{
    auto min = std::numeric_limits<double>::max();
    double umin = 0;
    for (size_t i = 0; i < c_.cpoints().size(); ++i) {
        auto u = itarg(i);
        if (auto cur = f(u); cur < min) {
            min = cur;
            umin = u;
        }
    }
    return {umin, min};
}

bool DistanceCurve::is_candidate(double d) const noexcept
{
    return std::any_of(std::begin(c_.cpoints()), std::end(c_.cpoints()),
                       [&d](auto& wp) {
                           auto p = pget(wp);
                           return cmp::le(p, d);
                       });
}

std::vector<SurfPoint> DistanceCurve::point_hull(double d) const
{
    auto s = c_.cpoints().size();

    std::vector<SurfPoint> values;
    values.reserve(s);

    for (decltype(s) i = 0; i < s; ++i) {
        auto& p = c_.cpoints()[i];
        auto u = itarg(i);
        auto v = pget(p) - d;
        values.emplace_back(u, v);
    }

    return graham_scan(values);
}

bool DistanceCurve::eliminate_segment(double d) noexcept
{
    static constexpr auto npos = size_t(-1);

    auto convex_hull = point_hull(d);
    auto roots = single_eliminate(convex_hull, c_.pfront(), c_.pback());

    for (size_t i = roots.size() - 1; i != npos; --i) {
        if (auto v = tocparg(roots[i], bool(i % 2)); f(v) < d) {
            auto it = std::begin(roots) + (i + 1);
            roots.erase(it);
        }
    }

    auto result = !roots.empty();
    if (result) {
        c_.refine_knots(roots);
        auto pc = c_.bezier_patches();

        if (roots.size() == 1
            && cmp::le(convex_hull.front().u, roots.front())) {
            c_ = pc.at(0);
        } else {
            c_ = pc.at(1);
        }
    }

    return result;
}

bool DistanceCurve::peak_point() const noexcept
{
    bool flag = false;
    bool result = true;
    auto s = c_.cpoints().size();
    double prev, cur = pget(c_.cpoints().front());

    for (size_t i = 1; i < s; ++i) {
        prev = cur;
        cur = pget(c_.cpoints()[i]);

        if (prev > cur) {
            if (flag) {
                result = false;
                break;
            } else {
                result = i;
            }
        } else if (!cmp::near(prev, cur)) {
            if (!flag) {
                result = i;
                flag = true;
            }
        } else {
            result = false;
            break;
        }
    }

    return result;
}

} // namespace gm
