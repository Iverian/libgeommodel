#include <gm/oriented_edge.hpp>

#include <fmt/ostream.h>

#include <util/to_string.hpp>

using namespace std;

namespace gm {

OrientedEdge::OrientedEdge(const Edge& edge, bool orientation)
    : edge_(edge)
    , orient_(orientation)
{
}

OrientedEdge::OrientedEdge(Edge&& edge, bool orientation) noexcept
    : edge_(edge)
    , orient_(orientation)
{
}

const Edge& OrientedEdge::OrientedEdge::edge() const noexcept
{
    return edge_;
}

bool OrientedEdge::OrientedEdge::orientation() const noexcept
{
    return orient_;
}

ostream& operator<<(ostream& os, const OrientedEdge& oedge)
{
    fmt::print(os, "{{ \"edge\": {0}, \"orientation\": {1} }}", oedge.edge_,
               oedge.orient_);
    return os;
}

} // namespace gm