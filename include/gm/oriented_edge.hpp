#ifndef GEOM_MODEL_INCLUDE_GM_ORIENTED_EDGE_HPP_
#define GEOM_MODEL_INCLUDE_GM_ORIENTED_EDGE_HPP_

#include "edge.hpp"
#include "exports.hpp"

namespace gm {

class GM_EXPORT OrientedEdge {
public:
    OrientedEdge(const Edge& edge, bool orientation);
    OrientedEdge(Edge&& edge, bool orientation) noexcept;

    [[nodiscard]] const Edge& edge() const noexcept;
    [[nodiscard]] bool orientation() const noexcept;

    friend std::ostream& operator<<(std::ostream& os,
                                    const OrientedEdge& oedge);

private:
    Edge edge_;
    bool orient_;
};

} // namespace gm

#endif // GEOM_MODEL_INCLUDE_GM_ORIENTED_EDGE_HPP_
