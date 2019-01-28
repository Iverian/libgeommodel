#ifndef GEOM_MODEL_INCLUDE_GM_BSPLINE_CURVE_H_
#define GEOM_MODEL_INCLUDE_GM_BSPLINE_CURVE_H_

#include "abstract_curve.h"

#include <vector>

namespace gm {

class BSplineCurve : public AbstractCurve {
public:
    struct Impl;

    BSplineCurve(size_t degree, std::vector<double> knots,
                 std::vector<Point> points,
                 std::vector<double> weights = std::vector<double>());
    BSplineCurve(size_t degree, const std::vector<size_t>& knot_multiplies,
                 const std::vector<double>& knot_list,
                 std::vector<Point> points,
                 std::vector<double> weights = std::vector<double>());

    Point f(double u) const noexcept override;
    Vec df(double u) const noexcept override;
    Vec df2(double u) const noexcept override;
    double project(const Point& p) const override;

    double pfront() const noexcept;
    double pback() const noexcept;

protected:
    std::ostream& print(std::ostream& os) const override;

private:
    std::shared_ptr<Impl> pimpl_;
};

} // namespace gm

#endif // GEOM_MODEL_INCLUDE_GM_BSPLINE_CURVE_H_
