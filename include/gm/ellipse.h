#ifndef GEOM_MODEL_INCLUDE_ELLIPSE_H_
#define GEOM_MODEL_INCLUDE_ELLIPSE_H_

#include "abstract_curve.h"
#include "axis.h"

namespace gm {

class Ellipse : public AbstractCurve {
public:
    Ellipse() noexcept;
    Ellipse(double rx, double ry, Axis ax) noexcept;

    double rx() const noexcept;
    double ry() const noexcept;
    const Axis& ax() const noexcept;

    Point f(double u) const noexcept override;
    Vec df(double u) const noexcept override;
    Vec df2(double u) const noexcept override;
    double project(const Point& p) const override;
    std::optional<double> project_greater(const Point& p, double min) const
        noexcept override;

protected:
    std::ostream& print(std::ostream& os) const override;

private:
    double rx_;
    double ry_;
    Axis ax_;
};

} // namespace gm

#endif // GEOM_MODEL_INCLUDE_ELLIPSE_H_
