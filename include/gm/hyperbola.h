#ifndef GEOM_MODEL_INCLUDE_GM_HYPERBOLA_H_
#define GEOM_MODEL_INCLUDE_GM_HYPERBOLA_H_

#include "abstract_curve.h"
#include "axis.h"

namespace gm {

class Hyperbola : public AbstractCurve {
public:
    Hyperbola() noexcept;
    Hyperbola(double rx, double ry, Axis ax) noexcept;

    Point f(double u) const noexcept override;
    Vec df(double u) const noexcept override;
    Vec df2(double u) const noexcept override;
    double project(const Point& p) const override;

protected:
    std::ostream& print(std::ostream& os) const override;

private:
    double rx_;
    double ry_;
    Axis ax_;
};

} // namespace gm

#endif // GEOM_MODEL_INCLUDE_GM_HYPERBOLA_H_
