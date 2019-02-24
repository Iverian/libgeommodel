#ifndef GEOM_MODEL_SRC_GEOM_BSPLINE_CURVE_IMPL_H_
#define GEOM_MODEL_SRC_GEOM_BSPLINE_CURVE_IMPL_H_

#include <bspline/basic_bspline_curve.h>
#include <gm/bspline_curve.h>

#include <memory>
#include <optional>
#include <ostream>
#include <vector>

class BSplineCurveProjector;

namespace gm {
struct CurveProjector;

class BSplineCurve::Impl {
public:
    static constexpr size_t dim = 3;
    using BezierPatch = ::BezierPatch<dim>;
    using Super = ::BasicBsplineCurve<dim>;

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

    double pfront() const noexcept;
    double pback() const noexcept;

    std::vector<BezierPatch> bezier_patches() const;

private:
    void init_curve(size_t order, const std::vector<double>& k,
                    const std::vector<Point>& p, const std::vector<double>& w);

    Super c_;
    mutable std::unique_ptr<CurveProjector> proj_;
};

} // namespace gm

#endif // GEOM_MODEL_SRC_GEOM_BSPLINE_CURVE_IMPL_H_