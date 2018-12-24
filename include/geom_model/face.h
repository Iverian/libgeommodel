#ifndef GEOM_MODEL_INCLUDE_FACE_H_
#define GEOM_MODEL_INCLUDE_FACE_H_

#include "abstract_surface.h"
#include "oriented_edge.h"

#include <memory>
#include <vector>

class Face;

std::ostream& operator<<(std::ostream& os, const Face& f);

using FaceBound = std::vector<OrientedEdge>;

class Face {
public:
    Face(std::shared_ptr<AbstractSurface> surface = nullptr,
         bool same_sense = true);
    Face(std::shared_ptr<AbstractSurface> surface, bool same_sense,
         const FaceBound& outer_loop,
         const std::vector<FaceBound>& inner_loop = std::vector<FaceBound>());

    bool same_sense() const;
    const std::vector<FaceBound>& boundaries() const;
    const FaceBound& outer() const;
    const AbstractSurface& surface() const;

private:
    struct Impl;
    std::shared_ptr<Impl> pimpl_;
};

#endif // GEOM_MODEL_INCLUDE_FACE_H_
