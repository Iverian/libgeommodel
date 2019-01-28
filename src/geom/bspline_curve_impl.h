#ifndef GEOM_MODEL_SRC_GEOM_BSPLINE_CURVE_IMPL_H_
#define GEOM_MODEL_SRC_GEOM_BSPLINE_CURVE_IMPL_H_

#include <gm/bspline_curve.h>
#include <primitive/wpoint.h>

#include <memory>
#include <optional>
#include <ostream>
#include <vector>

#include "cox_de_boor.h"

class BSplineCurveProjector;

namespace gm {

struct CurveProjector;

struct BezierPatch {
    explicit BezierPatch(size_t order = 0);
    size_t order() const noexcept;

    double front;
    double back;
    std::vector<CPoint> cp;
};

struct BSplineCurve::Impl {

    ~Impl();
    Impl(Impl&&) noexcept;
    Impl& operator=(Impl&&) noexcept;
    Impl(const Impl& other);
    Impl& operator=(const Impl& other);

    explicit Impl(const BezierPatch& patch);
    Impl(size_t degree, std::vector<double> knots, std::vector<Point> points,
         std::vector<double> weights = std::vector<double>());
    Impl(size_t degree, const std::vector<size_t>& knot_mult,
         const std::vector<double>& knot_list, std::vector<Point> points,
         std::vector<double> weights = std::vector<double>());

    Point f(double u) const noexcept;
    Vec df(double u) const noexcept;
    Vec df2(double u) const noexcept;
    std::ostream& print(std::ostream& os) const;

    double project(const Point& p) const;
    double ustep(double u) const;
    std::optional<double> newton_iter(const Point& p, double init,
                                      size_t max_iter) const;
    bool init_between(const Point& p, double a, double b) const;
    double bound_check(double u) const;

    double pfront() const noexcept;
    double pback() const noexcept;

    std::vector<BezierPatch> bezier_patches() const;

private:
    void init_cpoints(const std::vector<Point>& p,
                      const std::vector<double>& w);
    void init_cdb();

    size_t order_;
    std::vector<double> knots_;
    std::vector<CPoint> cpoints_;
    CoxDeBoor<CPoint> cdb_;
    mutable std::unique_ptr<CurveProjector> proj_;
};

} // namespace gm

#endif // GEOM_MODEL_SRC_GEOM_BSPLINE_CURVE_IMPL_H_