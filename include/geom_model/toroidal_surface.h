#ifndef GEOM_MODEL_INCLUDE_TOROIDAL_SURFACE_H
#define GEOM_MODEL_INCLUDE_TOROIDAL_SURFACE_H

#include "abstract_surface.h"
#include "axis.h"

#include <memory>

class ToroidalSurface : public AbstractSurface {
public:
    ~ToroidalSurface() override;
    ToroidalSurface(ToroidalSurface&&) noexcept;
    ToroidalSurface& operator=(ToroidalSurface&&) noexcept;
    ToroidalSurface(const ToroidalSurface& other);
    ToroidalSurface& operator=(const ToroidalSurface& other);

    ToroidalSurface();
    ToroidalSurface(double r0, double r1, Axis ax);

    Point f(const ParametricPoint& p) const override;
    Vec dfu(const ParametricPoint& p) const override;
    Vec dfv(const ParametricPoint& p) const override;
    Vec dfuu(const ParametricPoint& p) const override;
    Vec dfvv(const ParametricPoint& p) const override;
    Vec dfuv(const ParametricPoint& p) const override;
    std::ostream& print(std::ostream& os) const override;
    ParametricPoint project(const Point& p) const override;

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl_;
};

#endif // GEOM_MODEL_INCLUDE_TOROIDAL_SURFACE_H
