#ifndef GEOM_MODEL_SRC_GEOM_CURVE_PROJECTOR_H_
#define GEOM_MODEL_SRC_GEOM_CURVE_PROJECTOR_H_

#include <bspline/bspline_curve_impl.h>
#include <bspline/distance_curve.h>
#include <util/debug.h>

#include <vector>

namespace gm {

struct CurveProjector {
public:
    explicit CurveProjector(const BSplineCurve::Impl::Super& impl);
    double call(const Point& p) const;

    std::optional<double>
    minimize(const BSplineCurve::Impl::Super::CPoint::Proj& p, double u0,
             const double& pfront, const double& pback) const noexcept;
    double bord_check(double u, const double& pfront,
                      const double& pback) const noexcept;
    double armijo_step(double u, double h, const Vec& w, const Vec& d,
                       const BSplineCurve::Impl::Super::CPoint::Proj& p) const
        noexcept;

private:
    const BSplineCurve::Impl::Super* parent_;
    std::vector<BSplineCurve::Impl::BezierPatch> patches_;
};

} // namespace gm

#endif // GEOM_MODEL_SRC_GEOM_CURVE_PROJECTOR_H_