#include <gm/cylindrical_surface.hpp>

#include <util/math.hpp>

#include <fmt/ostream.h>

namespace gm {

CylindricalSurface::CylindricalSurface() noexcept
    : r_(1)
    , ax_()
{
}

CylindricalSurface::CylindricalSurface(double r, Axis ax) noexcept
    : r_(r)
    , ax_(std::move(ax))
{
}

Point CylindricalSurface::f(const SurfPoint& p) const noexcept
{
    return ax_.pglobal(r_ * ::cos(p.u), r_ * ::sin(p.u), r_ * p.v);
}

Vec CylindricalSurface::dfu(const SurfPoint& p) const noexcept
{
    return r_ * ax_.vglobal(-::sin(p.u), ::cos(p.u), 0);
}

Vec CylindricalSurface::dfv(const SurfPoint& p) const noexcept
{
    return ax_.vglobal(0, 0, r_);
}

Vec CylindricalSurface::dfuu(const SurfPoint& p) const noexcept
{
    return -r_ * ax_.vglobal(::cos(p.u), ::sin(p.u), 0);
}

Vec CylindricalSurface::dfuv(const SurfPoint& p) const noexcept
{
    return ax_.vglobal(0, 0, 0);
}

Vec CylindricalSurface::dfvv(const SurfPoint& p) const noexcept
{
    return ax_.vglobal(0, 0, 0);
}

std::ostream& CylindricalSurface::print(std::ostream& os) const
{
    fmt::print(os,
               "{{ \"type\": \"cylindrical\", \"r\": {0}, \"axis\": {1} }}",
               r_, ax_);
    return os;
}

SurfPoint CylindricalSurface::project(const Point& p) const
{
    auto [c, x, y, z] = ax_.get_view();
    auto v = dot(p - c, z) / r_;
    auto w = (p - c) - z * dot(p - c, z);
    auto [first, second] = atan2v(dot(w, y), dot(w, x));

    return (dist(f({first, v}), p) < dist(f({second, v}), p))
        ? SurfPoint {first, v}
        : SurfPoint {second, v};
}

} // namespace gm