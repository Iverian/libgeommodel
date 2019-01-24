#include <gm/conical_surface.h>
#include <gm/line.h>

#include <util/math.h>

#include <fmt/ostream.h>

using namespace std;

namespace gm {

ConicalSurface::ConicalSurface() noexcept
    : r_(1)
    , ta_(tan(M_PI / 6))
    , ax_()
{
}

ConicalSurface::ConicalSurface(double r, double a, Axis ax) noexcept
    : r_(r)
    , ta_(tan(a))
    , ax_(move(ax))
{
}

Point ConicalSurface::f(const SurfPoint& p) const noexcept
{
    return ax_.pglobal((r_ + p.v * ta_) * ::cos(p.u),
                       (r_ + p.v * ta_) * ::sin(p.u), p.v);
}

Vec ConicalSurface::dfu(const SurfPoint& p) const noexcept
{
    return ax_.vglobal(-(r_ + p.v * ta_) * ::sin(p.u),
                       (r_ + p.v * ta_) * ::cos(p.u), 0);
}

Vec ConicalSurface::dfv(const SurfPoint& p) const noexcept
{
    return ax_.vglobal(ta_ * ::cos(p.u), ta_ * ::sin(p.u), 1);
}

Vec ConicalSurface::dfuu(const SurfPoint& p) const noexcept
{
    return ax_.vglobal(-(r_ + p.v * ta_) * ::cos(p.u),
                       (r_ + p.v * ta_) * ::sin(p.u), 0);
}

Vec ConicalSurface::dfuv(const SurfPoint& p) const noexcept
{
    return ax_.vglobal(-ta_ * ::sin(p.u), ta_ * ::cos(p.u), 0);
}

Vec ConicalSurface::dfvv(const SurfPoint& p) const noexcept
{
    return ax_.vglobal(0, 0, 0);
}

ostream& ConicalSurface::print(ostream& os) const
{
    fmt::print(os,
               "{{ \"type\": \"conical_surface\", \"r\": {0}, \"ta\": {1}, "
               "\"axis\": {2} }}",
               r_, ta_, ax_);
    return os;
}

SurfPoint ConicalSurface::project(const Point& p) const
{
    auto [c, x, y, z] = ax_.get_view();
    auto w = (p - c) - z * dot(p - c, z);
    auto u = atan2v(dot(w, y), dot(w, x));
    array<double, 2> v;

    transform(cbegin(u), cend(u), begin(v), [&](const auto& i) {
        return Line(ta_ * ::cos(i) * x + ta_ * ::sin(i) * y + z,
                    r_ * ::cos(i) * x + r_ * ::sin(i) * y + c)
            .project(p);
    });

    return dist(f({u[0], v[0]}), p) < dist(f({u[1], v[1]}), p)
        ? SurfPoint {u[0], v[0]}
        : SurfPoint {u[1], v[1]};
}

} // namespace gm