#ifndef GEOM_MODEL_INCLUDE_TOROIDAL_SURFACE_H
#define GEOM_MODEL_INCLUDE_TOROIDAL_SURFACE_H

#include "abstract_surface.h"
#include "axis.h"

#include <memory>

namespace gm {

class ToroidalSurface : public AbstractSurface {
public:
    ToroidalSurface() noexcept;
    ToroidalSurface(double r0, double r1, Axis ax) noexcept;

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
    double r0_;
    double r1_;
    Axis ax_;
};

} // namespace gm

#endif // GEOM_MODEL_INCLUDE_TOROIDAL_SURFACE_H
