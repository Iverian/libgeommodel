#include <gm/axis.hpp>
#include <gm/ellipse.hpp>

#include <util/math.hpp>

#include <fmt/ostream.h>

namespace gm {

Ellipse::Ellipse() noexcept
    : rx_(1)
    , ry_(1)
    , ax_()
{
}

Ellipse::Ellipse(double rx, double ry, Axis ax) noexcept
    : rx_(rx)
    , ry_(ry)
    , ax_(std::move(ax))
{
}

double Ellipse::rx() const noexcept
{
    return rx_;
}

double Ellipse::ry() const noexcept
{
    return ry_;
}

const Axis& Ellipse::ax() const noexcept
{
    return ax_;
}

Point Ellipse::f(double u) const noexcept
{
    return ax_.pglobal(rx_ * ::cos(u), ry_ * ::sin(u), 0);
}

Vec Ellipse::df(double u) const noexcept
{
    return ax_.vglobal(-rx_ * ::sin(u), ry_ * ::cos(u), 0);
}

Vec Ellipse::df2(double u) const noexcept
{
    return ax_.vglobal(-rx_ * ::cos(u), -ry_ * ::sin(u), 0);
}

double Ellipse::project(const Point& p) const
{
    auto [c, x, y, z] = ax_.get_view();
    auto w = (p - c) - z * dot(p - c, z);
    auto [first, second] = atan2v(rx_ * dot(w, y), ry_ * dot(w, x));

    return dist(f(first), p) < dist(f(second), p) ? first : second;
}

std::optional<double> Ellipse::project_greater(const Point& x,
                                               double min) const noexcept
{
    auto result = project(x);
    while (result < min) {
        result += 2 * M_PI;
    }
    return result;
}

std::ostream& Ellipse::print(std::ostream& os) const
{
    fmt::print(
        os,
        "{{ \"type\": \"ellipse\", \"rx\": {0}, \"ry\": {1}, \"axis\": {2} }}",
        rx_, ry_, ax_);
    return os;
}

} // namespace gm