#ifndef GEOM_MODEL_INCLUDE_BSPLINE_SURFACE_H_
#define GEOM_MODEL_INCLUDE_BSPLINE_SURFACE_H_

#include "abstract_surface.h"

#include <memory>
#include <vector>

namespace gm {

class BSplineSurface : public AbstractSurface {
public:
    struct Impl;

    ~BSplineSurface() override;
    BSplineSurface(BSplineSurface&&) noexcept;
    BSplineSurface& operator=(BSplineSurface&&) noexcept;
    BSplineSurface(const BSplineSurface& other);
    BSplineSurface& operator=(const BSplineSurface& other);

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

    Point f(const SurfPoint& p) const noexcept override;
    Vec dfu(const SurfPoint& p) const noexcept override;
    Vec dfv(const SurfPoint& p) const noexcept override;
    Vec dfuu(const SurfPoint& p) const noexcept override;
    Vec dfuv(const SurfPoint& p) const noexcept override;
    Vec dfvv(const SurfPoint& p) const noexcept override;

    SurfPoint project(const Point& p) const override;

protected:
    std::ostream& print(std::ostream& os) const override;

private:
    std::unique_ptr<Impl> pimpl_;
};

} // namespace gm

#endif // GEOM_MODEL_INCLUDE_BSPLINE_SURFACE_H_
