#ifndef GEOM_MODEL_INCLUDE_GM_HYPERBOLA_HPP_
#define GEOM_MODEL_INCLUDE_GM_HYPERBOLA_HPP_

#include "abstract_curve.h"
#include "axis.h"
#include "exports.h"

namespace gm {

class GM_EXPORT Hyperbola : public AbstractCurve {
public:
    Hyperbola() noexcept;
    Hyperbola(double rx, double ry, Axis ax) noexcept;

    [[nodiscard]] Point f(double u) const noexcept override;
    [[nodiscard]] Vec df(double u) const noexcept override;
    [[nodiscard]] Vec df2(double u) const noexcept override;
    [[nodiscard]] double project(const Point& p) const override;

protected:
    std::ostream& print(std::ostream& os) const override;

private:
    double rx_;
    double ry_;
    Axis ax_;
};

} // namespace gm

#endif // GEOM_MODEL_INCLUDE_GM_HYPERBOLA_HPP_
