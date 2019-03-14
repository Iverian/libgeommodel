#ifndef GEOM_MODEL_SRC_GEOM_DISTANCE_CURVE_H_
#define GEOM_MODEL_SRC_GEOM_DISTANCE_CURVE_H_

#include <bspline/basic_bspline_curve.h>
#include <bspline/bspline_curve_impl.h>
#include <gm/surf_point.h>

namespace gm {

class DistanceCurve {
public:
    using Super = ::BasicBsplineCurve<1>;

    DistanceCurve();
    DistanceCurve(const BSplineCurve::Impl::BezierPatch& patch,
                  const Point& p);

    double pfront() const noexcept;
    double pback() const noexcept;
    double itarg(size_t cp_index) const noexcept;
    double argti(double arg) const noexcept;
    double tocparg(double arg, bool dir = false) const noexcept;

    double f(double u) const noexcept;

    std::pair<double, double> min_init() const noexcept;
    bool is_candidate(double d) const noexcept;
    std::vector<SurfPoint> point_hull(double d) const;
    bool eliminate_segment(double d) noexcept;
    bool peak_point() const noexcept;

private:
    Super c_;
};

} // namespace gm

#endif // GEOM_MODEL_SRC_GEOM_DISTANCE_CURVE_H_