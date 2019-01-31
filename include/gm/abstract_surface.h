#ifndef GEOM_MODEL_INCLUDE_GM_ABSTRACT_SURFACE_H_
#define GEOM_MODEL_INCLUDE_GM_ABSTRACT_SURFACE_H_

#include "point.h"
#include "surf_point.h"
#include "vec.h"

#include <functional>
#include <iostream>
#include <memory>
#include <optional>

namespace gm {
class Plane;

class AbstractSurface : public std::enable_shared_from_this<AbstractSurface> {
public:
    virtual ~AbstractSurface();
    [[nodiscard]] virtual Point f(const SurfPoint& p) const noexcept = 0;
    [[nodiscard]] virtual SurfPoint project(const Point& p) const = 0;

    [[nodiscard]] virtual Vec dfu(const SurfPoint& p) const noexcept;
    [[nodiscard]] virtual Vec dfv(const SurfPoint& p) const noexcept;
    [[nodiscard]] virtual Vec dfuu(const SurfPoint& p) const noexcept;
    [[nodiscard]] virtual Vec dfvv(const SurfPoint& p) const noexcept;
    [[nodiscard]] virtual Vec dfuv(const SurfPoint& p) const noexcept;

    [[nodiscard]] Point operator()(const SurfPoint& p) const noexcept;
    [[nodiscard]] Vec normal(const SurfPoint& p) const noexcept;
    [[nodiscard]] Vec unit_normal(const SurfPoint& p) const noexcept;
    [[nodiscard]] Plane tangent(const SurfPoint& p) const noexcept;
    [[nodiscard]] Point gproject(const Point& p) const noexcept;

    [[nodiscard]] std::function<Vec(double)> u_fixed(double v) const;
    [[nodiscard]] std::function<Vec(double)> v_fixed(double u) const;

    friend std::ostream& operator<<(std::ostream& os,
                                    const AbstractSurface& s);

protected:
    virtual std::ostream& print(std::ostream& os) const = 0;
};

} // namespace gm

#endif // GEOM_MODEL_INCLUDE_GM_ABSTRACT_SURFACE_H_
