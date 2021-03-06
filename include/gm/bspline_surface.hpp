#ifndef GEOM_MODEL_INCLUDE_GM_BSPLINE_SURFACE_HPP_
#define GEOM_MODEL_INCLUDE_GM_BSPLINE_SURFACE_HPP_

#include "abstract_surface.hpp"
#include "exports.hpp"

#include <memory>
#include <vector>

namespace gm {

class GM_EXPORT BSplineSurface : public AbstractSurface {
public:
    class Impl;

    BSplineSurface();
    BSplineSurface(size_t du, size_t dv, const std::vector<double>& ku,
                   const std::vector<double>& kv,
                   const std::vector<std::vector<Point>>& p,
                   const std::vector<std::vector<double>>& w
                   = std::vector<std::vector<double>>());

    BSplineSurface(size_t du, size_t dv, const std::vector<size_t>& ku_mult,
                   const std::vector<double>& ku_vals,
                   const std::vector<size_t>& kv_mult,
                   const std::vector<double>& kv_vals,
                   const std::vector<std::vector<Point>>& p,
                   const std::vector<std::vector<double>>& w
                   = std::vector<std::vector<double>>());

    [[nodiscard]] Point f(const SurfPoint& p) const noexcept override;
    [[nodiscard]] Vec dfu(const SurfPoint& p) const noexcept override;
    [[nodiscard]] Vec dfv(const SurfPoint& p) const noexcept override;
    [[nodiscard]] Vec dfuu(const SurfPoint& p) const noexcept override;
    [[nodiscard]] Vec dfuv(const SurfPoint& p) const noexcept override;
    [[nodiscard]] Vec dfvv(const SurfPoint& p) const noexcept override;

    [[nodiscard]] SurfPoint project(const Point& p) const override;

protected:
    std::ostream& print(std::ostream& os) const override;

private:
    std::shared_ptr<Impl> pimpl_;
};

} // namespace gm

#endif // GEOM_MODEL_INCLUDE_GM_BSPLINE_SURFACE_HPP_
