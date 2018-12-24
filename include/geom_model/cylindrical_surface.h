#ifndef GEOM_MODEL_INCLUDE_CYLINDRICAL_SURFACE_H_
#define GEOM_MODEL_INCLUDE_CYLINDRICAL_SURFACE_H_

#include "abstract_surface.h"
#include "axis.h"

#include <memory>

class CylindricalSurface : public AbstractSurface {
public:
    ~CylindricalSurface() override;
    CylindricalSurface(CylindricalSurface&&) noexcept;
    CylindricalSurface& operator=(CylindricalSurface&&) noexcept;
    CylindricalSurface(const CylindricalSurface& other);
    CylindricalSurface& operator=(const CylindricalSurface& other);

    CylindricalSurface();
    CylindricalSurface(double r, Axis ax);

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

#endif // GEOM_MODEL_INCLUDE_CYLINDRICAL_SURFACE_H_
