#include <gm/surf_point.hpp>

#include <fmt/ostream.h>

namespace gm {

SurfPoint abs(const SurfPoint& p)
{
    return {::abs(p.u), ::abs(p.v)};
}

std::ostream& operator<<(std::ostream& os, const SurfPoint& p)
{
    fmt::print(os, "[{0}, {1}]", p.u, p.v);
    return os;
}

} // namespace gm
