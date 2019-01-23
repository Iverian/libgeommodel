#include <gm/hyperbola.h>

#include <util/math.h>

#include <fmt/ostream.h>

using namespace std;

namespace gm {

Hyperbola::Hyperbola() noexcept
    : rx_(1)
    , ry_(1)
    , ax_()
{
}

Hyperbola::Hyperbola(double rx, double ry, Axis ax) noexcept
    : rx_(rx)
    , ry_(ry)
    , ax_(move(ax))
{
}

Point Hyperbola::f(double u) const noexcept
{
    return ax_.pglobal(rx_ * cosh(u), ry_ * sinh(u), 0);
}

Vec Hyperbola::df(double u) const noexcept
{
    return ax_.vglobal(rx_ * sinh(u), ry_ * cosh(u), 0);
}

Vec Hyperbola::df2(double u) const noexcept
{
    return ax_.vglobal(rx_ * cosh(u), ry_ * sinh(u), 0);
}

ostream& Hyperbola::print(ostream& os) const
{
    fmt::print(os,
               "{{ \"type\": \"hyperbola\", \"rx\": {0}, \"ry\": {1}, "
               "\"axis\": {2} }}",
               rx_, ry_, ax_);
    return os;
}

double Hyperbola::project(const Point& p) const
{
    auto [c, x, y, z] = ax_.get_view();
    auto w = (p - c) - z * dot(p - c, z);
    return atanh((rx_ * dot(w, y)) / (ry_ * dot(w, x)));
}

} // namespace gm