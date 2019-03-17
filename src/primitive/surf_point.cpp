#include <gm/surf_point.hpp>

#include <fmt/ostream.hpp>

using namespace std;

namespace gm {

SurfPoint abs(const SurfPoint& p)
{
    return {::abs(p.u), ::abs(p.v)};
}

ostream& operator<<(ostream& os, const SurfPoint& p)
{
    fmt::print(os, "[{0}, {1}]", p.u, p.v);
    return os;
}

} // namespace gm
