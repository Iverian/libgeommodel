#ifndef GEOM_MODEL_SRC_BSPLINE_DISTANCE_SURFACE_H_
#define GEOM_MODEL_SRC_BSPLINE_DISTANCE_SURFACE_H_

#include <bspline/basic_bspline_surface.h>
#include <bspline/bspline_surface_impl.h>
#include <gm/bspline_surface.h>
#include <gm/point.h>
#include <gm/surf_point.h>

#include <utility>
#include <vector>

namespace gm {

class DistanceSurface {
public:
    using Super = ::BasicBSplineSurface<1>;

    DistanceSurface();
    DistanceSurface(const BSplineSurface::Impl::BezierPatch& patch,
                    const Point& p);

    SurfPoint pfront() const noexcept;
    SurfPoint pback() const noexcept;

    SurfPoint itarg(const std::pair<size_t, size_t>& i) const noexcept;
    SurfPoint argti(SurfPoint arg) const noexcept;
    SurfPoint tocparg(SurfPoint arg,
                      std::pair<bool, bool> dir = {false, false}) const
        noexcept;

    double f(const SurfPoint& p) const noexcept;
    std::pair<SurfPoint, double> min_init() const noexcept;
    bool is_candidate(double d) const noexcept;
    std::vector<Point> point_hull(double d) noexcept;
    bool eliminate_segment(double d) noexcept;
    bool is_flat_enough() const noexcept;

private:
    Super c_;
};

} // namespace gm

#endif // GEOM_MODEL_SRC_BSPLINE_DISTANCE_SURFACE_H_