#ifndef GEOM_MODEL_INCLUDE_SPHERICAL_SURFACE_H_
#define GEOM_MODEL_INCLUDE_SPHERICAL_SURFACE_H_

#include "abstract_surface.h"
#include "axis.h"

#include <memory>

class SphericalSurface : public AbstractSurface {
public:
    ~SphericalSurface() override;
    SphericalSurface(SphericalSurface&&) noexcept;
    SphericalSurface& operator=(SphericalSurface&&) noexcept;
    SphericalSurface(const SphericalSurface& other);
    SphericalSurface& operator=(const SphericalSurface& other);

    SphericalSurface();
    SphericalSurface(double r, Axis ax);
    Vec dfu(const ParametricPoint& p) const override;
    Vec dfv(const ParametricPoint& p) const override;
    Vec dfuu(const ParametricPoint& p) const override;
    Vec dfvv(const ParametricPoint& p) const override;
    Vec dfuv(const ParametricPoint& p) const override;
    Point f(const ParametricPoint& p) const override;
    std::ostream& print(std::ostream& os) const override;
    ParametricPoint project(const Point& p) const override;

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl_;
};

#endif // GEOM_MODEL_INCLUDE_SPHERICAL_SURFACE_H_
