#include <gm/circle.hpp>
#include <gm/spherical_surface.hpp>

#include <util/math.hpp>

#include <fmt/ostream.hpp>

using namespace std;

namespace gm {

SphericalSurface::SphericalSurface() noexcept
    : r_(1)
    , ax_()
{
}

SphericalSurface::SphericalSurface(double r, Axis ax) noexcept
    : r_(r)
    , ax_(move(ax))
{
}

Point SphericalSurface::f(const SurfPoint& p) const noexcept
{
    return ax_.pglobal(r_ * ::cos(p.v) * ::cos(p.u), r_ * ::cos(p.v) * ::sin(p.u),
                       r_ * ::sin(p.v));
}

Vec SphericalSurface::dfu(const SurfPoint& p) const noexcept
{
    return ax_.vglobal(-r_ * ::cos(p.v) * ::sin(p.u), r_ * ::cos(p.v) * ::cos(p.u), 0);
}

Vec SphericalSurface::dfv(const SurfPoint& p) const noexcept
{
    return ax_.vglobal(-r_ * ::sin(p.v) * ::cos(p.u), -r_ * ::sin(p.v) * ::sin(p.u),
                       r_ * ::cos(p.v));
}

Vec SphericalSurface::dfuu(const SurfPoint& p) const noexcept
{
    return ax_.vglobal(-r_ * ::cos(p.v) * ::cos(p.u), -r_ * ::cos(p.v) * ::sin(p.u),
                       0);
}

Vec SphericalSurface::dfuv(const SurfPoint& p) const noexcept
{
    return ax_.vglobal(r_ * ::sin(p.v) * ::sin(p.u), -r_ * ::sin(p.v) * ::cos(p.u), 0);
}

Vec SphericalSurface::dfvv(const SurfPoint& p) const noexcept
{
    return ax_.vglobal(-r_ * ::cos(p.v) * ::cos(p.u), -r_ * ::cos(p.v) * ::sin(p.u),
                       -r_ * ::sin(p.v));
}

ostream& SphericalSurface::print(ostream& os) const
{
    fmt::print(os, "{{ \"type\": \"spherical\", \"r\": {0}, \"axis\": {1} }}",
               r_, ax_);
    return os;
}

SurfPoint SphericalSurface::project(const Point& p) const
{
    auto [c, x, y, z] = ax_.get_view();
    auto u = Circle(r_, ax_).project(p);
    auto v
        = Circle(r_, Axis::from_xy(x * ::cos(u) + y * ::sin(u), z, c)).project(p);
    return {u, v};
}

} // namespace gm
