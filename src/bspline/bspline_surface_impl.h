#ifndef GEOM_MODEL_SRC_GEOM_BSPLINE_SURFACE_IMPL_H_
#define GEOM_MODEL_SRC_GEOM_BSPLINE_SURFACE_IMPL_H_

#include <bspline/basic_bspline_surface.h>
#include <bspline/wpoint.h>
#include <gm/bspline_surface.h>

#include <array>
#include <memory>
#include <vector>

namespace gm {

class SurfaceProjector;

class BSplineSurface::Impl {
public:
    static constexpr size_t dim = 3;
    using BezierPatch = BezierSurfacePatch<dim>;
    using Super = BasicBSplineSurface<dim>;
    using OrderType = std::pair<size_t, size_t>;
    using KnotsType = Super::KnotsType;
    using CPointsType = Super::CPointsType;

    Impl();
    ~Impl();
    Impl(Impl&&) noexcept;
    Impl& operator=(Impl&&) noexcept;
    Impl(const Impl&);
    Impl& operator=(const Impl&);

    Impl(size_t du, size_t dv, const std::vector<double>& ku,
         const std::vector<double>& kv,
         const std::vector<std::vector<Point>>& p,
         const std::vector<std::vector<double>>& w
         = std::vector<std::vector<double>>());

    Impl(size_t du, size_t dv, const std::vector<size_t>& ku_mult,
         const std::vector<double>& ku_vals,
         const std::vector<size_t>& kv_mult,
         const std::vector<double>& kv_vals,
         const std::vector<std::vector<Point>>& p,
         const std::vector<std::vector<double>>& w
         = std::vector<std::vector<double>>());

    Point f(const SurfPoint& p) const noexcept;
    Vec dfu(const SurfPoint& p) const noexcept;
    Vec dfv(const SurfPoint& p) const noexcept;
    Vec dfuu(const SurfPoint& p) const noexcept;
    Vec dfvv(const SurfPoint& p) const noexcept;
    Vec dfuv(const SurfPoint& p) const noexcept;

    std::ostream& print(std::ostream& os) const;
    SurfPoint project(const Point& p) const;

    const OrderType& order() const noexcept;
    const KnotsType& knots() const noexcept;
    const CPointsType& cpoints() const noexcept;

private:
    void
    init_surface(std::pair<size_t, size_t> order,
                 const std::pair<std::vector<double>, std::vector<double>>& k,
                 const std::vector<std::vector<Point>>& p,
                 const std::vector<std::vector<double>>& w);

    void init_surface(std::pair<size_t, size_t> order,
                      std::pair<std::vector<double>, std::vector<double>>&& k,
                      const std::vector<std::vector<Point>>& p,
                      const std::vector<std::vector<double>>& w);

    Super c_;
    mutable std::unique_ptr<SurfaceProjector> proj_;
};

} // namespace gm

#endif // GEOM_MODEL_SRC_GEOM_BSPLINE_SURFACE_IMPL_H_