#ifndef GEOM_MODEL_SRC_GEOM_CURVE_PROJECTOR_H_
#define GEOM_MODEL_SRC_GEOM_CURVE_PROJECTOR_H_

#include <bspline/bspline_curve_impl.h>
#include <bspline/distance_curve.h>
#include <util/debug.h>

#include <vector>

namespace gm {

class CurveProjector {
public:
    explicit CurveProjector(const BSplineCurve::Impl& impl);
    double call(const Point& p) const;

    std::optional<double> minimize(const Point& p, double u0, const double& a,
                                   const double& b) const noexcept;
    double bord_check(double u, const double& a, const double& b) const
        noexcept;

private:
    const BSplineCurve::Impl* impl_;
    std::vector<BSplineCurve::Impl::BezierPatch> patches_;
};

} // namespace gm

#endif // GEOM_MODEL_SRC_GEOM_CURVE_PROJECTOR_H_