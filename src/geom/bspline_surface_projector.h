#ifndef GEOM_MODEL_SRC_GEOM_BSPLINE_SURFACE_PROJECTOR_H_
#define GEOM_MODEL_SRC_GEOM_BSPLINE_SURFACE_PROJECTOR_H_

#include <gm/point.h>
#include <gm/surf_point.h>

#include <primitive/wpoint.h>

class BSplineSurfaceProjector {
public:
    gm::SurfPoint operator()(const gm::Point& p) const;

private:
    size_t order_;

};

#endif // GEOM_MODEL_SRC_GEOM_BSPLINE_SURFACE_PROJECTOR_H_