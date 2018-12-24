#ifndef GEOM_MODEL_INCLUDE_ABSTRACT_CURVE_H_
#define GEOM_MODEL_INCLUDE_ABSTRACT_CURVE_H_

#include "point.h"
#include "vec.h"

#include <array>
#include <memory>
#include <optional>
#include <vector>

class AbstractCurve;

class AbstractCurve : public std::enable_shared_from_this<AbstractCurve> {
public:
    virtual ~AbstractCurve() = default;
    virtual Point f(double u) const = 0;
    virtual double project(const Point& p) const = 0;

    virtual Vec df(double u) const;
    virtual Vec df2(double u) const;
    virtual double length(double begin, double end) const;
    virtual std::optional<double> project_greater(const Point& p,
                                                  double min) const;

    Vec tangent(double u) const;
    Vec normal(double u) const;
    double curvature(double u) const;
    Point operator()(double u) const;
    Point gproject(const Point& p) const;

    friend std::ostream& operator<<(std::ostream& os, const AbstractCurve& c);

protected:
    virtual std::ostream& print(std::ostream& os) const = 0;

    std::optional<double> project_iterative(const Point& p, double init,
                                            size_t max_iter = 10000) const;
    bool is_init_in_interval(const Point& p,
                             const std::array<double, 2>& interval) const;
    double approx_length(double begin, double end, size_t mesh_size) const;
};

#endif // GEOM_MODEL_INCLUDE_ABSTRACT_CURVE_H_
