#include <gm/surf_point.hpp>

#include <fmt/ostream.h>

#include <cmath>
#include <utility>

namespace gm {

SurfPoint abs(const SurfPoint& p)
{
    return {std::abs(p.u), std::abs(p.v)};
}

double hypot(const SurfPoint& p)
{
    return std::hypot(p.u, p.v);
}

std::ostream& operator<<(std::ostream& os, const SurfPoint& p)
{
    fmt::print(os, "[{0}, {1}]", p.u, p.v);
    return os;
}

} // namespace gm
