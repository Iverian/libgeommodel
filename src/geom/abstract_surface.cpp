#include <gm/abstract_surface.h>
#include <gm/compare.h>
#include <gm/plane.h>
#include <gm/point.h>
#include <gm/surf_point.h>
#include <gm/vec.h>

#include <util/debug.h>
#include <util/math.h>

#include <algorithm>

using namespace std;

namespace gm {

Point AbstractSurface::operator()(const SurfPoint& p) const noexcept
{
    return f(p);
}

Vec AbstractSurface::normal(const SurfPoint& p) const noexcept
{
    return cross(dfu(p), dfv(p));
}

Vec AbstractSurface::unit_normal(const SurfPoint& p) const noexcept
{
    return unit(normal(p));
}

Plane AbstractSurface::tangent(const SurfPoint& p) const noexcept
{
    return Plane(Axis::from_xy(dfu(p), dfv(p), f(p)));
}

Point AbstractSurface::gproject(const Point& p) const noexcept
{
    return f(project(p));
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

function<Vec(double)> AbstractSurface::u_fixed(double v) const
{
    return [v, this](double u) { return Vec(f({u, v})); };
}

function<Vec(double)> AbstractSurface::v_fixed(double u) const
{
    return [u, this](double v) { return Vec(f({u, v})); };
}

ostream& operator<<(ostream& os, const AbstractSurface& s)
{
    return s.print(os);
}

} // namespace gm