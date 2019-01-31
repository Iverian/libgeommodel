#ifndef GEOM_MODEL_SRC_GEOM_BSPLINE_SURFACE_IMPL_H_
#define GEOM_MODEL_SRC_GEOM_BSPLINE_SURFACE_IMPL_H_

#include <gm/bspline_surface.h>
#include <gm/point.h>
#include <gm/surf_point.h>
#include <gm/vec.h>
#include <primitive/wpoint.h>

#include <array>
#include <vector>

#include "cox_de_boor.h"

namespace gm {

struct BSplineSurface::Impl {
public:
    using order_type = std::array<size_t, 2>;
    using knot_type = std::array<std::vector<double>, 2>;
    using cpoint_type = std::vector<CPoint>;
    using cox_de_boor_type = std::vector<CoxDeBoor<CPoint>>;

    struct CPointSize {
        size_t n;
        size_t m;
    };

    Impl();
    Impl(const Impl&) = default;
    Impl(Impl&&) noexcept = default;
    Impl& operator=(const Impl&) = default;
    Impl& operator=(Impl&&) noexcept = default;

    Impl(const order_type& order, const knot_type& knots,
         const cpoint_type& cpoints, const CPointSize& size);

    Impl(const order_type& order, const knot_type& knots,
         const std::vector<std::vector<CPoint>>& cpoints);

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

    const order_type& order() const noexcept;
    const knot_type& knots() const noexcept;
    const cpoint_type& cpoints() const noexcept;
    const CPointSize& size() const noexcept;

    Impl bezier_patches() const;

private:
    void init_cpoints(const std::vector<std::vector<CPoint>>& cp);
    void init_cpoints(const std::vector<std::vector<double>>& w,
                      const std::vector<std::vector<Point>>& p);
    void init_cdb();

    order_type order_;
    knot_type knots_;
    cpoint_type cpoints_;
    CPointSize size_;

    cox_de_boor_type cdb_;
};

} // namespace gm

#endif // GEOM_MODEL_SRC_GEOM_BSPLINE_SURFACE_IMPL_H_