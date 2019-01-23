#ifndef GEOM_MODEL_INCLUDE_PLANE_H_
#define GEOM_MODEL_INCLUDE_PLANE_H_

#include "abstract_surface.h"
#include "axis.h"

#include <memory>

namespace gm {

class Plane : public AbstractSurface {
public:
    Plane() noexcept;
    explicit Plane(Axis ax) noexcept;

    const Axis& ax() const noexcept;

    Point f(const SurfPoint& p) const noexcept override;
    Vec dfu(const SurfPoint& p) const noexcept override;
    Vec dfv(const SurfPoint& p) const noexcept override;
    Vec dfuu(const SurfPoint& p) const noexcept override;
    Vec dfvv(const SurfPoint& p) const noexcept override;
    Vec dfuv(const SurfPoint& p) const noexcept override;
    SurfPoint project(const Point& p) const override;

protected:
    std::ostream& print(std::ostream& os) const override;

private:
    Axis ax_;
};

} // namespace gm

#endif // GEOM_MODEL_INCLUDE_PLANE_H_
