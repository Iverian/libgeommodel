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

enum class WhereMin { NIL, FRONT, BACK };

class DistanceSurface {
public:
    using Super = ::BasicBSplineSurface<1>;

    DistanceSurface();
    DistanceSurface(const BSplineSurface::Impl::BezierPatch& patch,
                    const Point& p);

    SurfPoint pfront() const noexcept;
    SurfPoint pback() const noexcept;

    SurfPoint itarg(const std::pair<size_t, size_t>& i) const noexcept;
    SurfPoint itarg(size_t i) const noexcept;
    SurfPoint argti(SurfPoint arg) const noexcept;
    SurfPoint tocparg(SurfPoint arg,
                      std::pair<bool, bool> dir = {false, false}) const
        noexcept;

    double f(const SurfPoint& p) const noexcept;
    std::pair<SurfPoint, double> min_init() const noexcept;
    bool is_candidate(double d) const noexcept;
    std::pair<std::vector<gm::SurfPoint>, std::vector<gm::SurfPoint>>
    point_hull(double d) const noexcept;
    bool eliminate_segment(double d) noexcept;
    bool is_flat_enough() const noexcept;

    std::pair<WhereMin, WhereMin> min_bord_check() const noexcept;

private:
    Super c_;
};

} // namespace gm

#endif // GEOM_MODEL_SRC_BSPLINE_DISTANCE_SURFACE_HPP_