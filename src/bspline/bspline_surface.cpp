#include <bspline/bspline_surface_impl.hpp>
#include <gm/bspline_surface.hpp>

namespace gm {

BSplineSurface::BSplineSurface()
    : pimpl_(std::make_unique<Impl>())
{
}

BSplineSurface::BSplineSurface(size_t du, size_t dv,
                               const std::vector<double>& ku,
                               const std::vector<double>& kv,
                               const std::vector<std::vector<Point>>& p,
                               const std::vector<std::vector<double>>& w)
    : pimpl_(std::make_unique<Impl>(du, dv, ku, kv, p, w))
{
}

BSplineSurface::BSplineSurface(size_t du, size_t dv,
                               const std::vector<size_t>& ku_mult,
                               const std::vector<double>& ku_vals,
                               const std::vector<size_t>& kv_mult,
                               const std::vector<double>& kv_vals,
                               const std::vector<std::vector<Point>>& p,
                               const std::vector<std::vector<double>>& w)
    : pimpl_(std::make_unique<Impl>(du, dv, ku_mult, ku_vals, kv_mult, kv_vals,
                                    p, w))
{
}

Point BSplineSurface::f(const SurfPoint& p) const noexcept
{
    return pimpl_->f(p);
}

Vec BSplineSurface::dfu(const SurfPoint& p) const noexcept
{
    return pimpl_->dfu(p);
}

Vec BSplineSurface::dfv(const SurfPoint& p) const noexcept
{
    return pimpl_->dfv(p);
}

Vec BSplineSurface::dfuu(const SurfPoint& p) const noexcept
{
    return pimpl_->dfuu(p);
}

Vec BSplineSurface::dfuv(const SurfPoint& p) const noexcept
{
    return pimpl_->dfuv(p);
}

Vec BSplineSurface::dfvv(const SurfPoint& p) const noexcept
{
    return pimpl_->dfvv(p);
}

std::ostream& BSplineSurface::print(std::ostream& os) const
{
    return pimpl_->print(os);
}

SurfPoint BSplineSurface::project(const Point& p) const
{
    return pimpl_->project(p);
}

} // namespace gm
