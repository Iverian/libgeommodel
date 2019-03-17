#ifndef GEOM_MODEL_INCLUDE_GM_PARABOLA_HPP_
#define GEOM_MODEL_INCLUDE_GM_PARABOLA_HPP_

#include "abstract_curve.hpp"
#include "axis.hpp"
#include "exports.hpp"

namespace gm {

class GM_EXPORT Parabola : public AbstractCurve {
public:
    Parabola() noexcept;
    Parabola(double f, Axis ax) noexcept;

    [[nodiscard]] Point f(double u) const noexcept override;
    [[nodiscard]] Vec df(double u) const noexcept override;
    [[nodiscard]] Vec df2(double u) const noexcept override;
    [[nodiscard]] double project(const Point& p) const override;
    [[nodiscard]] double approx_length(double begin, double end,
                                       size_t n) const override;

protected:
    std::ostream& print(std::ostream& os) const override;

private:
    double f_;
    Axis ax_;
};

} // namespace gm

#endif // GEOM_MODEL_INCLUDE_GM_PARABOLA_HPP_
