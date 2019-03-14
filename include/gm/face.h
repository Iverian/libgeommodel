#ifndef GEOM_MODEL_INCLUDE_GM_FACE_H_
#define GEOM_MODEL_INCLUDE_GM_FACE_H_

#include "abstract_surface.h"
#include "exports.h"
#include "oriented_edge.h"

#include <memory>
#include <vector>

namespace gm {

using FaceBound = std::vector<OrientedEdge>;

class GM_EXPORT Face {
public:
    Face(std::shared_ptr<AbstractSurface> surface = nullptr,
         bool same_sense = true);
    Face(std::shared_ptr<AbstractSurface> surface, bool same_sense,
         const FaceBound& outer_loop,
         const std::vector<FaceBound>& inner_loop = std::vector<FaceBound>());

    [[nodiscard]] bool same_sense() const;
    [[nodiscard]] const std::vector<FaceBound>& boundaries() const;
    [[nodiscard]] const FaceBound& outer() const;
    [[nodiscard]] const AbstractSurface& surface() const;

private:
    struct Impl;
    std::shared_ptr<Impl> pimpl_;
};

GM_EXPORT std::ostream& operator<<(std::ostream& os, const Face& f);

} // namespace gm

#endif // GEOM_MODEL_INCLUDE_GM_FACE_H_
