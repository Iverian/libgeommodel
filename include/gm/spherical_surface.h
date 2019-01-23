#ifndef GEOM_MODEL_INCLUDE_SPHERICAL_SURFACE_H_
#define GEOM_MODEL_INCLUDE_SPHERICAL_SURFACE_H_

#include "abstract_surface.h"
#include "axis.h"

namespace gm {

class SphericalSurface : public AbstractSurface {
public:
    SphericalSurface() noexcept;
    SphericalSurface(double r, Axis ax) noexcept;

    Vec dfu(const SurfPoint& p) const noexcept override;
    Vec dfv(const SurfPoint& p) const noexcept override;
    Vec dfuu(const SurfPoint& p) const noexcept override;
    Vec dfvv(const SurfPoint& p) const noexcept override;
    Vec dfuv(const SurfPoint& p) const noexcept override;
    Point f(const SurfPoint& p) const noexcept override;
    SurfPoint project(const Point& p) const override;

protected:
    std::ostream& print(std::ostream& os) const override;

private:
    double r_;
    Axis ax_;
};

} // namespace gm

#endif // GEOM_MODEL_INCLUDE_SPHERICAL_SURFACE_H_
