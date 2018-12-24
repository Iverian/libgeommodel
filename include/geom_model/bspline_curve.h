#ifndef GEOM_MODEL_INCLUDE_BSPLINE_CURVE_H_
#define GEOM_MODEL_INCLUDE_BSPLINE_CURVE_H_

#include "abstract_curve.h"

#include <vector>

class BSplineCurve : public AbstractCurve {
public:
    ~BSplineCurve() override;
    BSplineCurve(BSplineCurve&&) noexcept;
    BSplineCurve& operator=(BSplineCurve&&) noexcept;
    BSplineCurve(const BSplineCurve& other);
    BSplineCurve& operator=(const BSplineCurve& other);

    BSplineCurve(size_t degree, std::vector<double> knots,
                 std::vector<Point> control_points,
                 std::vector<double> weights = std::vector<double>());
    BSplineCurve(size_t degree, const std::vector<size_t>& knot_multiplies,
                 const std::vector<double>& knot_list,
                 std::vector<Point> control_points,
                 std::vector<double> weights = std::vector<double>());

    Point f(double u) const override;
    Vec df(double u) const override;
    Vec df2(double u) const override;
    double project(const Point& p) const override;

    const double& param_front() const;
    const double& param_back() const;

    size_t get_order() const;
    const std::vector<double>& get_knots() const;
    const std::vector<Point>& get_control_points() const;
    const std::vector<double>& get_weights() const;

protected:
    std::ostream& print(std::ostream& os) const override;

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl_;
};
#endif // GEOM_MODEL_INCLUDE_BSPLINE_CURVE_H_
