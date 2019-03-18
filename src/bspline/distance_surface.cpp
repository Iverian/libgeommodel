#include "C:\Users\trololo\Documents\Projects\libgeommodel\include\gm\compare.hpp"
#include <bspline/distance_surface.hpp>
#include <bspline/util.hpp>
#include <bspline/wpoint.hpp>
#include <gm/surf_point.hpp>
#include <stdexcept>
#include <util/cyclic_iterator.hpp>
#include <util/math.hpp>

#include <algorithm>
#include <limits>
#include <utility>
#include <xutility>

#define pget(cp) ((cp).p()[0])

namespace gm {

DistanceSurface::DistanceSurface()
    : c_()
{
}

DistanceSurface::DistanceSurface(
    const BSplineSurface::Impl::BezierPatch& patch, const Point& p)
    : c_()
{
    auto order = std::make_pair(2 * patch.order().first - 1,
                                2 * patch.order().second - 1);
    auto deg
        = std::make_pair(patch.order().first - 1, patch.order().second - 1);
    auto knots = make_pair(std::vector<double>(2 * order.first),
                           std::vector<double>(2 * order.second));
    Super::CPointsType cp(order.first * order.second);
    auto r = p.raw();

    size_t i, j;
    double c;
    for (j = 0; j < order.second; ++j) {
        knots.second[j] = patch.pfront().v;
        knots.second[order.second + j] = patch.pback().v;
    }
    for (i = 0; i < order.first; ++i) {
        knots.first[i] = patch.pfront().u;
        knots.first[order.first + i] = patch.pback().u;

        auto i_first = (i < deg.first) ? size_t(0) : (i - deg.first);
        auto i_last = std::min(deg.first, i) + 1;
        for (j = 0; j < order.second; ++j) {
            auto& t = cp[j + order.second * i];

            auto j_first = (j < deg.second) ? size_t(0) : (j - deg.second);
            auto j_last = std::min(deg.second, j) + 1;

            for (auto p = i_first; p < i_last; ++p) {
                for (auto q = j_first; q < j_last; ++q) {
                    auto a = patch[{p, q}] - r;
                    auto b = patch[{i - p, j - q}] - r;

                    c = binom(deg.first, p) * binom(deg.first, i - p)
                        * binom(deg.second, q) * binom(deg.second, j - q);
                    t += c * wdot(a, b);
                }
            }
            c = binom(2 * deg.first, i) * binom(2 * deg.second, j);
            t /= c;
        }
    }
}

SurfPoint DistanceSurface::pfront() const noexcept
{
    return c_.pfront();
}

SurfPoint DistanceSurface::pback() const noexcept
{
    return c_.pback();
}

SurfPoint DistanceSurface::itarg(const std::pair<size_t, size_t>& i) const
    noexcept
{
    auto u = pfront().u
        + i.first * (pback().u - pfront().u) / (c_.order().first - 1);
    auto v = pfront().v
        + i.second * (pback().v - pfront().v) / (c_.order().second - 1);

    return {u, v};
}

SurfPoint DistanceSurface::itarg(size_t i) const noexcept
{
    auto& s = c_.shape();
    auto p = i / s.second;
    auto q = i % s.second;
    return itarg({p, q});
}

SurfPoint DistanceSurface::argti(SurfPoint arg) const noexcept
{
    auto u = (c_.order().first - 1) * (arg.u - pfront().u)
        / (pback().u - pfront().u);
    auto v = (c_.order().second - 1) * (arg.v - pfront().v)
        / (pback().v - pfront().v);

    return {u, v};
}
SurfPoint DistanceSurface::tocparg(SurfPoint arg,
                                   std::pair<bool, bool> dir) const noexcept
{
    arg = argti(arg);
    return itarg({dir.first ? ceil<size_t>(arg.u) : floor<size_t>(arg.u),
                  dir.second ? ceil<size_t>(arg.v) : floor<size_t>(arg.v)});
}

double DistanceSurface::f(const SurfPoint& p) const noexcept
{
    return c_.f(p)[0];
}

std::pair<SurfPoint, double> DistanceSurface::min_init() const noexcept
{
    auto s = c_.shape();
    auto dmin = std::numeric_limits<double>::max();
    SurfPoint argmin {};

    for (size_t i = 0; i < s.first; ++i) {
        for (size_t j = 0; j < s.second; ++j) {
            auto r = itarg({i, j});
            auto d = f(r);
            if (d < dmin) {
                argmin = r;
                dmin = d;
            }
        }
    }

    return {argmin, dmin};
}
bool DistanceSurface::is_candidate(double d) const noexcept
{
    return std::any_of(std::begin(c_.cpoints()), std::end(c_.cpoints()),
                       [&d](auto& wp) {
                           auto p = pget(wp);
                           return cmp::le(p, d);
                       });
}

std::pair<std::vector<gm::SurfPoint>, std::vector<gm::SurfPoint>>
DistanceSurface::point_hull(double d) const noexcept
{
    std::pair<std::vector<gm::SurfPoint>, std::vector<gm::SurfPoint>> result;
    auto s = c_.cpoints().size();
    result.first.reserve(s);
    result.second.reserve(s);
    for (size_t i = 0; i < s; ++i) {
        auto t = itarg(i);
        auto f = pget(c_.cpoints()[i]) - d;
        result.first.emplace_back(t.u, f);
        result.second.emplace_back(t.v, f);
    }

    return {graham_scan(result.first), graham_scan(result.second)};
}

bool DistanceSurface::eliminate_segment(double d) noexcept
{
    auto [cu, cv] = point_hull(d);
    auto ru = single_eliminate(cu, c_.pfront().u, c_.pback().u);
    auto rv = single_eliminate(cv, c_.pfront().v, c_.pback().v);

    auto result = !ru.empty() && !rv.empty();
    if (result) {
        c_.refine_knots({ru, rv});
        auto p = c_.bezier_patches();
        size_t i = 1, j = 1;
        if (ru.size() == 1 && cmp::le(cu.front().u, ru.front())) {
            i = 0;
        }
        if (rv.size() == 1 && cmp::le(cv.front().u, rv.front())) {
            j = 0;
        }
        c_ = p[i * (rv.size() + 1) + j];
    }

    return result;
}

bool DistanceSurface::is_flat_enough() const noexcept
{
    throw std::runtime_error("not implemented");
    return false;
}

} // namespace gm

#undef pget