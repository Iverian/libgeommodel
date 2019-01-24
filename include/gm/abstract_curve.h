#ifndef GEOM_MODEL_INCLUDE_GM_ABSTRACT_CURVE_H_
#define GEOM_MODEL_INCLUDE_GM_ABSTRACT_CURVE_H_

#include "point.h"
#include "vec.h"

#include <array>
#include <memory>
#include <optional>
#include <vector>

namespace gm {

class AbstractCurve : public std::enable_shared_from_this<AbstractCurve> {
public:
    virtual ~AbstractCurve();
    virtual Point f(double u) const noexcept = 0;
    virtual double project(const Point& p) const = 0;

    virtual Vec df(double u) const noexcept;
    virtual Vec df2(double u) const noexcept;
    virtual std::optional<double> project_greater(const Point& p,
                                                  double min) const noexcept;
    virtual double approx_length(double begin, double end, size_t n) const;

    Vec tangent(double u) const noexcept;
    Vec normal(double u) const noexcept;
    double curvature(double u) const noexcept;
    Point operator()(double u) const noexcept;
    Point gproject(const Point& p) const noexcept;

    friend std::ostream& operator<<(std::ostream& os, const AbstractCurve& c);

protected:
    virtual std::ostream& print(std::ostream& os) const = 0;
};

} // namespace gm

#endif // GEOM_MODEL_INCLUDE_GM_ABSTRACT_CURVE_H_
