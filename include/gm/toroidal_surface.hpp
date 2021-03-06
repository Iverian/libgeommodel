#ifndef GEOM_MODEL_INCLUDE_TOROIDAL_SURFACE_H
#define GEOM_MODEL_INCLUDE_TOROIDAL_SURFACE_H

#include "abstract_surface.hpp"
#include "axis.hpp"
#include "exports.hpp"

#include <memory>

namespace gm {

class GM_EXPORT ToroidalSurface : public AbstractSurface {
public:
    ToroidalSurface() noexcept;
    ToroidalSurface(double r0, double r1, Axis ax) noexcept;

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
    double r0_;
    double r1_;
    Axis ax_;
};

} // namespace gm

#endif // GEOM_MODEL_INCLUDE_TOROIDAL_SURFACE_H
