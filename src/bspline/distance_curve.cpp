#include <bspline/distance_curve.h>
#include <bspline/wpoint.h>
#include <util/cyclic_iterator.h>
#include <util/math.h>

#include <iterator>
#include <optional>

using namespace std;

double binom(size_t n, size_t k) __GM_NOEXCEPT_RELEASE__;
::vector<gm::SurfPoint> graham_scan(::vector<gm::SurfPoint>& points);
double polar_angle(const gm::SurfPoint& lhs, const gm::SurfPoint& rhs);
bool counter_clockwise(const gm::SurfPoint& a, const gm::SurfPoint& b,
                       const gm::SurfPoint& c);
::optional<double>
zero_intersect(const ::pair<gm::SurfPoint, gm::SurfPoint>& line);

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
    vector<double> knots(2 * order);
    vector<Super::CPoint> cp(order);

    auto pr = point.raw();
    auto deg = patch.order() - 1;
    double c = 0;

    for (size_t i = 0; i < order; ++i) {
        knots[i] = patch.pfront();
        knots[order + i] = patch.pback();

        auto& q = cp[i];
        auto first = (i < deg) ? size_t(0) : (i - deg);
        auto last = ::min(deg, i) + 1;
        for (auto k = first; k < last; ++k) {
            auto a = patch.cpoints()[k] - pr;
            auto b = patch.cpoints()[i - k] - pr;

            c = binom(deg, k) * binom(deg, i - k);
            q += c * wdot(a, b);
        }
        c = binom(2 * deg, i);
        q /= c;
    }

    c_ = Super(order, move(knots), move(cp));
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

double DistanceCurve::df(double u) const noexcept
{
    return c_.df(u)[0];
}

double DistanceCurve::df2(double u) const noexcept
{
    return c_.df2(u)[0];
}

::pair<double, double> DistanceCurve::min_init() const noexcept
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
    return ::any_of(::begin(c_.cpoints()), ::end(c_.cpoints()),
                    [&d](auto& wp) {
                        auto p = pget(wp);
                        return cmp::le(p, d);
                    });
}

std::vector<SurfPoint> DistanceCurve::point_hull(double d) const
{
    auto s = c_.cpoints().size();

    ::vector<SurfPoint> values;
    values.reserve(s);

    for (decltype(s) i = 0; i < s; ++i) {
        auto& p = c_.cpoints()[i];
        auto u = itarg(i);
        auto v = pget(p) - d;
        values.emplace_back(u, v);
    }

    return ::graham_scan(values);
}

bool DistanceCurve::eliminate_segment(double d) noexcept
{
    static constexpr auto max_roots = size_t(2);
    static constexpr auto npos = size_t(-1);

    auto convex_hull = point_hull(d);

    vector<double> roots;
    roots.reserve(max_roots);

    auto it = CyclicIterator(::begin(convex_hull), ::end(convex_hull));
    do {
        if (auto u = ::zero_intersect({*it, *::next(it)}); u) {
            auto v = u.value();
            if (!cmp::near(v, c_.pfront()) && !cmp::near(v, c_.pback())
                && (roots.empty() || !cmp::near(v, roots.back()))) {
                roots.emplace_back(v);
            }
        }
    } while (++it, roots.size() != max_roots && it.iter() != it.first());
    ::sort(::begin(roots), ::end(roots));

    for (size_t i = roots.size() - 1; i != npos; --i) {
        if (auto v = tocparg(roots[i], bool(i % 2)); f(v) < d) {
            auto it = ::begin(roots) + (i + 1);
            roots.erase(it);
        }
    }

    auto result = !roots.empty();
    if (result) {
        auto new_c = c_.refine_knots(roots);
        auto pc = new_c.bezier_patches();

        if (roots.size() == 1
            && cmp::le(convex_hull.front().u, roots.front())) {
            c_ = pc.at(0);
        } else {
            c_ = pc.at(1);
        }
    }

    return result;
}

::optional<size_t> DistanceCurve::peak_point() const noexcept
{
    bool flag = false;
    ::optional<size_t> result = 0;
    auto s = c_.cpoints().size();
    double prev, cur = pget(c_.cpoints().front());

    for (size_t i = 1; i < s; ++i) {
        prev = cur;
        cur = pget(c_.cpoints()[i]);

        if (prev > cur) {
            if (flag) {
                result = nullopt;
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
            result = nullopt;
            break;
        }
    }

    return result;
}

} // namespace gm

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

#define next_to_top(stack) (*::prev(::end(stack), 2))
#define top(stack) ((stack).back())

::vector<gm::SurfPoint> graham_scan(::vector<gm::SurfPoint>& points)
{
    ::vector<gm::SurfPoint> result;
    result.reserve(points.size());

    {
        auto min = ::min_element(
            ::begin(points), ::end(points),
            [](auto& lhs, auto& rhs) { return lhs.v < rhs.v; });
        ::swap(*::begin(points), *min);
    }
    ::sort(::next(::begin(points)), ::end(points),
           [& p = points[0]](auto& lhs, auto& rhs) {
               return polar_angle(lhs, p) < polar_angle(rhs, p);
           });
    result.push_back(points[0]);
    result.push_back(points[1]);

    for (size_t i = 2; i < points.size(); ++i) {
        auto& q = points[i];
        while (result.size() > 1
               && !counter_clockwise(next_to_top(result), top(result), q)) {
            result.pop_back();
        }
        result.push_back(q);
    }

    result.shrink_to_fit();
    return result;
}

#undef next_to_top
#undef top

double polar_angle(const gm::SurfPoint& lhs, const gm::SurfPoint& rhs)
{
    auto p = rhs - lhs;
    return ::atan2(p.v, p.u);
}

bool counter_clockwise(const gm::SurfPoint& a, const gm::SurfPoint& b,
                       const gm::SurfPoint& c)
{
    auto value = (b.u - a.u) * (c.v - a.v) - (b.v - a.v) * (c.u - a.u);
    return value > 0;
}

::optional<double>
zero_intersect(const ::pair<gm::SurfPoint, gm::SurfPoint>& line)
{
    auto& [a, b] = line;

    if (gm::cmp::zero(a.v, gm::Tolerance::MAX)) {
        return a.u;
    }
    if (gm::cmp::zero(b.v, gm::Tolerance::MAX)) {
        return b.u;
    }
    if (::signbit(a.v) != ::signbit(b.v)) {
        return (a.u * b.v - b.u * a.v) / (b.v - a.v);
    }

    return nullopt;
}