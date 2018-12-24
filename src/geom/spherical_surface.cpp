#include <geom_model/circle.h>
#include <geom_model/spherical_surface.h>

#include <util/math.h>

#include <fmt/ostream.h>

using namespace std;

struct SphericalSurface::Impl {
    Impl();
    Impl(double r, Axis ax);

    Point f(const ParametricPoint& p) const;
    Vec dfu(const ParametricPoint& p) const;
    Vec dfv(const ParametricPoint& p) const;
    Vec dfuu(const ParametricPoint& p) const;
    Vec dfuv(const ParametricPoint& p) const;
    Vec dfvv(const ParametricPoint& p) const;
    ostream& print(ostream& os) const;
    ParametricPoint project(const Point& p) const;

private:
    double r_;
    Axis ax_;
};

SphericalSurface::~SphericalSurface() = default;
SphericalSurface::SphericalSurface(SphericalSurface&&) noexcept = default;
SphericalSurface& SphericalSurface::operator=(SphericalSurface&&) noexcept
    = default;

SphericalSurface::SphericalSurface(const SphericalSurface& other)
    : pimpl_(make_unique<Impl>(*other.pimpl_))
{
}

SphericalSurface& SphericalSurface::operator=(const SphericalSurface& other)
{
    *pimpl_ = *other.pimpl_;
    return *this;
}

SphericalSurface::SphericalSurface()
    : pimpl_(make_unique<Impl>())
{
}

SphericalSurface::SphericalSurface(double r, Axis ax)
    : pimpl_(make_unique<Impl>(r, move(ax)))
{
}

Point SphericalSurface::f(const ParametricPoint& p) const
{
    return pimpl_->f(p);
}

Vec SphericalSurface::dfu(const ParametricPoint& p) const
{
    return pimpl_->dfu(p);
}

Vec SphericalSurface::dfv(const ParametricPoint& p) const
{
    return pimpl_->dfv(p);
}

Vec SphericalSurface::dfuu(const ParametricPoint& p) const
{
    return pimpl_->dfuu(p);
}

Vec SphericalSurface::dfvv(const ParametricPoint& p) const
{
    return pimpl_->dfvv(p);
}

Vec SphericalSurface::dfuv(const ParametricPoint& p) const
{
    return pimpl_->dfuv(p);
}

ostream& SphericalSurface::print(ostream& os) const
{
    return pimpl_->print(os);
}

ParametricPoint SphericalSurface::project(const Point& p) const
{
    return pimpl_->project(p);
}

SphericalSurface::Impl::Impl()
    : r_(1)
    , ax_()
{
}

SphericalSurface::Impl::Impl(double r, Axis ax)
    : r_(r)
    , ax_(move(ax))
{
}

Point SphericalSurface::Impl::f(const ParametricPoint& p) const
{
    return ax_.pglobal(r_ * cos(p.v) * cos(p.u), r_ * cos(p.v) * sin(p.u),
                       r_ * sin(p.v));
}

Vec SphericalSurface::Impl::dfu(const ParametricPoint& p) const
{
    return ax_.vglobal(-r_ * cos(p.v) * sin(p.u), r_ * cos(p.v) * cos(p.u), 0);
}

Vec SphericalSurface::Impl::dfv(const ParametricPoint& p) const
{
    return ax_.vglobal(-r_ * sin(p.v) * cos(p.u), -r_ * sin(p.v) * sin(p.u),
                       r_ * cos(p.v));
}

Vec SphericalSurface::Impl::dfuu(const ParametricPoint& p) const
{
    return ax_.vglobal(-r_ * cos(p.v) * cos(p.u), -r_ * cos(p.v) * sin(p.u),
                       0);
}

Vec SphericalSurface::Impl::dfuv(const ParametricPoint& p) const
{
    return ax_.vglobal(r_ * sin(p.v) * sin(p.u), -r_ * sin(p.v) * cos(p.u), 0);
}

Vec SphericalSurface::Impl::dfvv(const ParametricPoint& p) const
{
    return ax_.vglobal(-r_ * cos(p.v) * cos(p.u), -r_ * cos(p.v) * sin(p.u),
                       -r_ * sin(p.v));
}

ostream& SphericalSurface::Impl::print(ostream& os) const
{
    fmt::print(os, "{{ \"type\": \"spherical\", \"r\": {0}, \"axis\": {1} }}",
               r_, ax_);
    return os;
}

ParametricPoint SphericalSurface::Impl::project(const Point& p) const
{
    auto [c, x, y, z] = ax_.get_view();
    auto u = Circle(r_, ax_).project(p);
    auto v
        = Circle(r_, Axis::from_xy(x * cos(u) + y * sin(u), z, c)).project(p);
    return {u, v};
}