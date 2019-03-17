#ifndef GEOM_MODEL_SRC_BSPLINE_SURFACE_PROJECTOR_HPP_
#define GEOM_MODEL_SRC_BSPLINE_SURFACE_PROJECTOR_HPP_

#include <bspline/bspline_surface_impl.h>
#include <gm/surf_point.h>
#include <util/debug.h>

#include <vector>

namespace gm {

class SurfaceProjector {
public:
    explicit SurfaceProjector(const BSplineSurface::Impl& impl);
    SurfPoint call(const Point& p) const;

    std::optional<SurfPoint> minimize(const Point& p, SurfPoint r,
                                      const SurfPoint& a,
                                      const SurfPoint& b) const noexcept;
    SurfPoint next_step(SurfPoint r, const Vec& w, const Vec& fu,
                        const Vec& fv) const noexcept;
    SurfPoint bord_check(SurfPoint r, const SurfPoint& a,
                         const SurfPoint& b) const noexcept;

private:
    const BSplineSurface::Impl* impl_;
};

} // namespace gm

#endif // GEOM_MODEL_SRC_BSPLINE_SURFACE_PROJECTOR_HPP_