#include <gm/circle.hpp>

#include <util/math.hpp>

#include <fmt/ostream.h>

namespace gm {

Circle::Circle() noexcept
    : Ellipse()
{
}

Circle::Circle(double r, Axis ax) noexcept
    : Ellipse(r, r, std::move(ax))
{
}

double Circle::r() const noexcept
{
    return rx();
}

double Circle::approx_length(double begin, double end, size_t n) const
{
    return r() * fabs(end - begin);
}

std::ostream& Circle::print(std::ostream& os) const
{
    fmt::print(os, "{{ \"type\": \"circle\", \"r\": {0}, \"axis\": {1} }}",
               r(), ax());
    return os;
}

} // namespace gm