#ifndef GEOM_MODEL_INCLUDE_GM_ORIENTED_EDGE_H_
#define GEOM_MODEL_INCLUDE_GM_ORIENTED_EDGE_H_

#include "edge.h"

namespace gm {

class OrientedEdge {
public:
    OrientedEdge(const Edge& edge, bool orientation);
    OrientedEdge(Edge&& edge, bool orientation) noexcept;

    const Edge& edge() const noexcept;
    bool orientation() const noexcept;

    friend std::ostream& operator<<(std::ostream& os,
                                    const OrientedEdge& oedge);

private:
    Edge edge_;
    bool orient_;
};

} // namespace gm

#endif // GEOM_MODEL_INCLUDE_GM_ORIENTED_EDGE_H_
