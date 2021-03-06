#ifndef GEOM_MODEL_INCLUDE_GM_PLANE_HPP_
#define GEOM_MODEL_INCLUDE_GM_PLANE_HPP_

#include "abstract_surface.hpp"
#include "axis.hpp"
#include "exports.hpp"

#include <memory>

namespace gm {

class GM_EXPORT Plane : public AbstractSurface {
public:
    Plane() noexcept;
    explicit Plane(Axis ax) noexcept;

    [[nodiscard]] const Axis& ax() const noexcept;

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
    Axis ax_;
};

} // namespace gm

#endif // GEOM_MODEL_INCLUDE_GM_PLANE_HPP_
