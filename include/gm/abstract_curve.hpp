#ifndef GEOM_MODEL_INCLUDE_GM_ABSTRACT_CURVE_HPP_
#define GEOM_MODEL_INCLUDE_GM_ABSTRACT_CURVE_HPP_

#include "exports.hpp"
#include "point.hpp"
#include "vec.hpp"

#include <array>
#include <memory>
#include <optional>
#include <vector>

namespace gm {

class GM_EXPORT AbstractCurve
    : public std::enable_shared_from_this<AbstractCurve> {
public:
    virtual ~AbstractCurve();
    [[nodiscard]] virtual Point f(double u) const noexcept = 0;
    [[nodiscard]] virtual double project(const Point& p) const = 0;

    [[nodiscard]] virtual Vec df(double u) const noexcept;
    [[nodiscard]] virtual Vec df2(double u) const noexcept;
    [[nodiscard]] virtual std::optional<double>
    project_greater(const Point& p, double min) const noexcept;
    [[nodiscard]] virtual double approx_length(double begin, double end,
                                               size_t n) const;

    [[nodiscard]] Vec tangent(double u) const noexcept;
    [[nodiscard]] Vec normal(double u) const noexcept;
    [[nodiscard]] Vec unit_normal(double u) const __GM_NOEXCEPT_RELEASE__;
    [[nodiscard]] double curvature(double u) const noexcept;
    [[nodiscard]] Point operator()(double u) const noexcept;
    [[nodiscard]] Point gproject(const Point& p) const noexcept;

    GM_EXPORT friend std::ostream& operator<<(std::ostream& os,
                                              const AbstractCurve& c);

protected:
    virtual std::ostream& print(std::ostream& os) const = 0;
};

} // namespace gm

#endif // GEOM_MODEL_INCLUDE_GM_ABSTRACT_CURVE_HPP_
