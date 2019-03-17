#include <gm/circle.hpp>
#include <gm/toroidal_surface.hpp>

#include <util/math.hpp>

#include <fmt/ostream.h>

namespace gm {

ToroidalSurface::ToroidalSurface() noexcept
    : r0_(0.5)
    , r1_(1)
    , ax_()
{
}

ToroidalSurface::ToroidalSurface(double r0, double r1, Axis ax) noexcept
    : r0_(r0)
    , r1_(r1)
    , ax_(std::move(ax))
{
}

Point ToroidalSurface::f(const SurfPoint& p) const noexcept
{
    return ax_.pglobal((r1_ + r0_ * ::cos(p.v)) * ::cos(p.u),
                       (r1_ + r0_ * ::cos(p.v)) * ::sin(p.u),
                       r0_ * ::sin(p.v));
}

Vec ToroidalSurface::dfu(const SurfPoint& p) const noexcept
{
    return (r1_ + r0_ * ::cos(p.v)) * ax_.vglobal(-::sin(p.u), ::cos(p.u), 0);
}

Vec ToroidalSurface::dfv(const SurfPoint& p) const noexcept
{
    return r0_
        * ax_.vglobal(::sin(p.v) * ::cos(p.u), ::sin(p.v) * ::sin(p.u),
                      ::cos(p.v));
}

Vec ToroidalSurface::dfuu(const SurfPoint& p) const noexcept
{
    return -(r1_ + r0_ * ::cos(p.v)) * ax_.vglobal(::cos(p.u), ::sin(p.u), 0);
}

Vec ToroidalSurface::dfuv(const SurfPoint& p) const noexcept
{
    return r0_ * ::sin(p.v) * ax_.vglobal(::sin(p.u), ::cos(p.u), 0);
}

Vec ToroidalSurface::dfvv(const SurfPoint& p) const noexcept
{
    return -r0_
        * ax_.vglobal(::cos(p.v) * ::cos(p.u), ::cos(p.v) * ::sin(p.u),
                      ::sin(p.v));
}

std::ostream& ToroidalSurface::print(std::ostream& os) const
{
    fmt::print(
        os,
        "{{ \"type\": \"toroidal\", \"r0\":{0}, \"r1\":{1}, \"axis\":{2} }}",
        r0_, r1_, ax_);
    return os;
}

SurfPoint ToroidalSurface::project(const Point& p) const
{
    auto [c, x, y, z] = ax_.get_view();
    auto u = Circle((r0_ + r1_) / 2, ax_).project(p);
    auto v = Circle(r0_,
                    Axis::from_xy(::cos(u) * x + ::sin(u) * y, z,
                                  c + r1_ * ::cos(u) * x + r1_ * ::sin(u) * y))
                 .project(p);
    return SurfPoint {u, v};
}

} // namespace gm