#ifndef GEOM_MODEL_INCLUDE_CONICAL_SURFACE_H_
#define GEOM_MODEL_INCLUDE_CONICAL_SURFACE_H_

#include "abstract_surface.h"
#include "axis.h"

#include <memory>

class ConicalSurface : public AbstractSurface {
public:
    ~ConicalSurface() override;
    ConicalSurface(ConicalSurface&&) noexcept;
    ConicalSurface& operator=(ConicalSurface&&) noexcept;
    ConicalSurface(const ConicalSurface& other);
    ConicalSurface& operator=(const ConicalSurface& other);

    ConicalSurface();
    ConicalSurface(double r, double a, Axis ax);

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

#endif // GEOM_MODEL_INCLUDE_CONICAL_SURFACE_H_
