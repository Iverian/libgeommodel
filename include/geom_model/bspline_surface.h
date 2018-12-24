#ifndef GEOM_MODEL_INCLUDE_BSPLINE_SURFACE_H_
#define GEOM_MODEL_INCLUDE_BSPLINE_SURFACE_H_

#include "abstract_surface.h"

#include <memory>
#include <vector>

class BSplineSurface : public AbstractSurface {
public:
    ~BSplineSurface() override;
    BSplineSurface(BSplineSurface&&) noexcept;
    BSplineSurface& operator=(BSplineSurface&&) noexcept;
    BSplineSurface(const BSplineSurface& other);
    BSplineSurface& operator=(const BSplineSurface& other);

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

    Point f(const ParametricPoint& p) const override;
    Vec dfu(const ParametricPoint& p) const override;
    Vec dfv(const ParametricPoint& p) const override;
    Vec dfuu(const ParametricPoint& p) const override;
    Vec dfuv(const ParametricPoint& p) const override;
    Vec dfvv(const ParametricPoint& p) const override;

    std::ostream& print(std::ostream& os) const override;
    ParametricPoint project(const Point& p) const override;

private:
    struct Impl;
    std::shared_ptr<Impl> pimpl_;
};

#endif // GEOM_MODEL_INCLUDE_BSPLINE_SURFACE_H_
