#include <gm/plane.hpp>

#include <fmt/ostream.h>

namespace gm {

Plane::Plane() noexcept
    : ax_()
{
}

Plane::Plane(Axis ax) noexcept
    : ax_(std::move(ax))
{
}

const Axis& Plane::ax() const noexcept
{
    return ax_;
}

Point Plane::f(const SurfPoint& p) const noexcept
{
    return ax_.pglobal(p.u, p.v, 0);
}

Vec Plane::dfu(const SurfPoint& p) const noexcept
{
    return ax_.vglobal(1, 0, 0);
}

Vec Plane::dfv(const SurfPoint& p) const noexcept
{
    return ax_.vglobal(0, 1, 0);
}

Vec Plane::dfuu(const SurfPoint& p) const noexcept
{
    return ax_.vglobal(0, 0, 0);
}

Vec Plane::dfuv(const SurfPoint& p) const noexcept
{
    return ax_.vglobal(0, 0, 0);
}

Vec Plane::dfvv(const SurfPoint& p) const noexcept
{
    return ax_.vglobal(0, 0, 0);
}

std::ostream& Plane::print(std::ostream& os) const
{
    fmt::print(os, "{{ \"type\": \"plane\", \"axis\": {0} }}", ax_);
    return os;
}

SurfPoint Plane::project(const Point& p) const
{
    auto [c, x, y, z] = ax_.get_view();
    auto w = (p - c) - z * dot(p - c, z);
    return {dot(w, x), dot(w, y)};
}

} // namespace gm