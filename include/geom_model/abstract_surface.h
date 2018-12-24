#ifndef GEOM_MODEL_INCLUDE_ABSTRACT_SURFACE_H_
#define GEOM_MODEL_INCLUDE_ABSTRACT_SURFACE_H_

#include "parametric_point.h"
#include "point.h"
#include "vec.h"

#include <functional>
#include <iostream>
#include <memory>
#include <optional>

class AbstractSurface;
class Plane;

std::ostream& operator<<(std::ostream& os, const AbstractSurface& s);

class AbstractSurface : public std::enable_shared_from_this<AbstractSurface> {
public:
    virtual ~AbstractSurface() = default;
    virtual Point f(const ParametricPoint& p) const = 0;
    virtual std::ostream& print(std::ostream& os) const = 0;
    virtual ParametricPoint project(const Point& p) const = 0;

    virtual Vec dfu(const ParametricPoint& p) const;
    virtual Vec dfv(const ParametricPoint& p) const;
    virtual Vec dfuu(const ParametricPoint& p) const;
    virtual Vec dfvv(const ParametricPoint& p) const;
    virtual Vec dfuv(const ParametricPoint& p) const;

    Point operator()(const ParametricPoint& p) const;
    Vec normal(const ParametricPoint& p) const;
    Plane tangent(const ParametricPoint& p) const;
    Point gproject(const Point& p) const;

    std::function<Vec(double)> u_fixed(double v) const;
    std::function<Vec(double)> v_fixed(double u) const;

protected:
    std::optional<ParametricPoint>
    project_iterative(const Point& p, const ParametricPoint& init,
                      size_t max_iter = 10000) const;
    bool is_init_in_square(const Point& p,
                           const std::array<ParametricPoint, 2>& diag) const;

private:
    ParametricPoint project_to_step(const Point& x,
                                    const ParametricPoint& p) const;
};

#endif // GEOM_MODEL_INCLUDE_ABSTRACT_SURFACE_H_
