#ifndef GEOM_MODEL_INCLUDE_GM_SPHERICAL_SURFACE_HPP_
#define GEOM_MODEL_INCLUDE_GM_SPHERICAL_SURFACE_HPP_

#include "abstract_surface.hpp"
#include "axis.hpp"
#include "exports.hpp"

namespace gm {

class GM_EXPORT SphericalSurface : public AbstractSurface {
public:
    SphericalSurface() noexcept;
    SphericalSurface(double r, Axis ax) noexcept;

    [[nodiscard]] Vec dfu(const SurfPoint& p) const noexcept override;
    [[nodiscard]] Vec dfv(const SurfPoint& p) const noexcept override;
    [[nodiscard]] Vec dfuu(const SurfPoint& p) const noexcept override;
    [[nodiscard]] Vec dfvv(const SurfPoint& p) const noexcept override;
    [[nodiscard]] Vec dfuv(const SurfPoint& p) const noexcept override;
    [[nodiscard]] Point f(const SurfPoint& p) const noexcept override;
    [[nodiscard]] SurfPoint project(const Point& p) const override;

protected:
    std::ostream& print(std::ostream& os) const override;

private:
    double r_;
    Axis ax_;
};

} // namespace gm

#endif // GEOM_MODEL_INCLUDE_GM_SPHERICAL_SURFACE_HPP_
