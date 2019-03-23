#ifndef GEOM_MODEL_SRC_BSPLINE_DISTANCE_SURFACE_HPP_
#define GEOM_MODEL_SRC_BSPLINE_DISTANCE_SURFACE_HPP_

#include <bspline/basic_bspline_surface.hpp>
#include <bspline/bspline_surface_impl.hpp>
#include <gm/bspline_surface.hpp>
#include <gm/point.hpp>
#include <gm/surf_point.hpp>

#include <utility>
#include <vector>

namespace gm {

class DistanceSurface {
public:
    using Super = ::BasicBSplineSurface<1>;

    DistanceSurface();
    DistanceSurface(const BSplineSurface::Impl::BezierPatch& patch,
                    const Point& p);

    const std::pair<size_t, size_t>& order() const noexcept;
    const std::pair<size_t, size_t>& shape() const noexcept;
    SurfPoint pfront() const noexcept;
    SurfPoint pback() const noexcept;

    SurfPoint itarg(const std::pair<size_t, size_t>& i) const noexcept;
    SurfPoint itarg(size_t i) const noexcept;

    bool is_candidate(double d) const noexcept;
    std::pair<std::vector<gm::SurfPoint>, std::vector<gm::SurfPoint>>
    point_hull(double d) const noexcept;
    bool eliminate_segment(double d) noexcept;

private:
    Super c_;
};

} // namespace gm

#endif // GEOM_MODEL_SRC_BSPLINE_DISTANCE_SURFACE_HPP_