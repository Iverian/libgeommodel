#include <gm/circle.h>

#include <util/math.h>

#include <fmt/ostream.h>

using namespace std;

namespace gm {

Circle::Circle() noexcept
    : Ellipse()
{
}

Circle::Circle(double r, Axis ax) noexcept
    : Ellipse(r, r, move(ax))
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

ostream& Circle::print(ostream& os) const
{
    fmt::print(os, "{{ \"type\": \"circle\", \"r\": {0}, \"axis\": {1} }}",
               r(), ax());
    return os;
}

} // namespace gm