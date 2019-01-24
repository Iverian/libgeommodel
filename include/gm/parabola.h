#ifndef GEOM_MODEL_INCLUDE_GM_PARABOLA_H_
#define GEOM_MODEL_INCLUDE_GM_PARABOLA_H_

#include "abstract_curve.h"
#include "axis.h"

namespace gm {

class Parabola : public AbstractCurve {
public:
    Parabola() noexcept;
    Parabola(double f, Axis ax) noexcept;

    Point f(double u) const noexcept override;
    Vec df(double u) const noexcept override;
    Vec df2(double u) const noexcept override;
    double project(const Point& p) const override;
    double approx_length(double begin, double end, size_t n) const override;

protected:
    std::ostream& print(std::ostream& os) const override;

private:
    double f_;
    Axis ax_;
};

} // namespace gm

#endif // GEOM_MODEL_INCLUDE_GM_PARABOLA_H_
