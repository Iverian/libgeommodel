#include <gm/abstract_surface.hpp>
#include <gm/compare.hpp>
#include <gm/plane.hpp>
#include <gm/point.hpp>
#include <gm/surf_point.hpp>
#include <gm/vec.hpp>

#include <util/math.hpp>

#include <algorithm>

namespace gm {

AbstractSurface::~AbstractSurface() = default;

Point AbstractSurface::operator()(const SurfPoint& p) const noexcept
{
    return f(p);
}

Vec AbstractSurface::normal(const SurfPoint& p) const noexcept
{
    return cross(dfu(p), dfv(p));
}

Vec AbstractSurface::unit_normal(const SurfPoint& p) const
    __GM_NOEXCEPT_RELEASE__
{
    return unit(normal(p));
}

Axis AbstractSurface::tangent(const SurfPoint& p) const noexcept
{
    return Axis::from_xy(dfu(p), dfv(p), f(p));
}

Point AbstractSurface::gproject(const Point& p) const noexcept
{
    return f(project(p));
}

[[nodiscard]] std::pair<double, double>
AbstractSurface::main_curvatures(const SurfPoint& p) const noexcept
{
    auto du = dfu(p);
    auto dv = dfv(p);
    auto duu = dfuu(p);
    auto duv = dfuv(p);
    auto dvv = dfuv(p);
    double g[] = {sqr(du), dot(du, dv), sqr(dv)};
    double gg = std::sqrt(g[0] * g[2] - sqr(g[1]));
    double b[] = {mixed(duu, du, dv) / gg, mixed(duv, du, dv) / gg,
                  mixed(dvv, du, dv) / gg};

    auto x = b[0] * g[2] - 2 * b[1] * g[1] + b[2] * g[0];
    auto y = std::sqrt(
        sqr(b[0] * g[2]) + sqr(b[2] * g[0]) - 2 * b[0] * b[2] * g[0] * g[2]
        + 4
            * (b[0] * b[2] * sqr(g[1]) + sqr(b[1]) * g[0] * g[2]
               - b[0] * b[1] * g[1] * g[2] - b[1] * b[2] * g[0] * g[1]));

    return {(x - y) / (2 * gg), (x + y) / (2 * gg)};
}

Vec AbstractSurface::dfu(const SurfPoint& p) const noexcept
{
    return ::diff<Vec>(u_fixed(p.v), p.u);
}

Vec AbstractSurface::dfv(const SurfPoint& p) const noexcept
{
    return ::diff<Vec>(v_fixed(p.u), p.v);
}

Vec AbstractSurface::dfuu(const SurfPoint& p) const noexcept
{
    return ::diff2<Vec>(u_fixed(p.v), p.u);
}

Vec AbstractSurface::dfvv(const SurfPoint& p) const noexcept
{
    return ::diff2<Vec>(v_fixed(p.u), p.v);
}

Vec AbstractSurface::dfuv(const SurfPoint& p) const noexcept
{
    return ::diff11<Vec>([this](SurfPoint t) { return Vec(f(t)); }, p);
}

std::function<Vec(double)> AbstractSurface::u_fixed(double v) const
{
    return [v, this](double u) { return Vec(f({u, v})); };
}

std::function<Vec(double)> AbstractSurface::v_fixed(double u) const
{
    return [u, this](double v) { return Vec(f({u, v})); };
}

std::ostream& operator<<(std::ostream& os, const AbstractSurface& s)
{
    return s.print(os);
}

} // namespace gm