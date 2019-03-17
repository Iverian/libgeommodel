#ifndef GEOM_MODEL_INCLUDE_GM_ELLIPSE_HPP_
#define GEOM_MODEL_INCLUDE_GM_ELLIPSE_HPP_

#include "abstract_curve.h"
#include "axis.h"
#include "exports.h"

namespace gm {

class GM_EXPORT Ellipse : public AbstractCurve {
public:
    Ellipse() noexcept;
    Ellipse(double rx, double ry, Axis ax) noexcept;

    [[nodiscard]] double rx() const noexcept;
    [[nodiscard]] double ry() const noexcept;
    [[nodiscard]] const Axis& ax() const noexcept;

    [[nodiscard]] Point f(double u) const noexcept override;
    [[nodiscard]] Vec df(double u) const noexcept override;
    [[nodiscard]] Vec df2(double u) const noexcept override;
    [[nodiscard]] double project(const Point& p) const override;
    [[nodiscard]] std::optional<double> project_greater(const Point& p,
                                                        double min) const
        noexcept override;

protected:
    std::ostream& print(std::ostream& os) const override;

private:
    double rx_;
    double ry_;
    Axis ax_;
};

} // namespace gm

#endif // GEOM_MODEL_INCLUDE_GM_ELLIPSE_HPP_
