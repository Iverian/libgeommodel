#ifndef GEOM_MODEL_INCLUDE_GM_CYLINDRICAL_SURFACE_H_
#define GEOM_MODEL_INCLUDE_GM_CYLINDRICAL_SURFACE_H_

#include "abstract_surface.h"
#include "axis.h"

namespace gm {

class CylindricalSurface : public AbstractSurface {
public:
    CylindricalSurface() noexcept;
    CylindricalSurface(double r, Axis ax) noexcept;

    [[nodiscard]] Point f(const SurfPoint& p) const noexcept override;
    [[nodiscard]] Vec dfu(const SurfPoint& p) const noexcept override;
    [[nodiscard]] Vec dfv(const SurfPoint& p) const noexcept override;
    [[nodiscard]] Vec dfuu(const SurfPoint& p) const noexcept override;
    [[nodiscard]] Vec dfvv(const SurfPoint& p) const noexcept override;
    [[nodiscard]] Vec dfuv(const SurfPoint& p) const noexcept override;

    [[nodiscard]] SurfPoint project(const Point& p) const override;

protected:
    std::ostream& print(std::ostream& os) const override;

private:
    double r_;
    Axis ax_;
};

} // namespace gm

#endif // GEOM_MODEL_INCLUDE_GM_CYLINDRICAL_SURFACE_H_
