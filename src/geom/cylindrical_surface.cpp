#include <geom_model/cylindrical_surface.h>

#include <util/math.h>

#include <fmt/ostream.h>

using namespace std;

struct CylindricalSurface::Impl {
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

CylindricalSurface::~CylindricalSurface() = default;
CylindricalSurface::CylindricalSurface(CylindricalSurface&&) noexcept
    = default;
CylindricalSurface& CylindricalSurface::
operator=(CylindricalSurface&&) noexcept
    = default;

CylindricalSurface::CylindricalSurface(const CylindricalSurface& other)
    : pimpl_(make_unique<Impl>(*other.pimpl_))
{
}

CylindricalSurface& CylindricalSurface::
operator=(const CylindricalSurface& other)
{
    *pimpl_ = *other.pimpl_;
    return *this;
}

CylindricalSurface::CylindricalSurface()
    : pimpl_(make_unique<Impl>())
{
}

CylindricalSurface::CylindricalSurface(double r, Axis ax)
    : pimpl_(make_unique<Impl>(r, ax))
{
}

Point CylindricalSurface::f(const ParametricPoint& p) const
{
    return pimpl_->f(p);
}

Vec CylindricalSurface::dfu(const ParametricPoint& p) const
{
    return pimpl_->dfu(p);
}

Vec CylindricalSurface::dfv(const ParametricPoint& p) const
{
    return pimpl_->dfv(p);
}

Vec CylindricalSurface::dfuu(const ParametricPoint& p) const
{
    return pimpl_->dfuu(p);
}

Vec CylindricalSurface::dfvv(const ParametricPoint& p) const
{
    return pimpl_->dfvv(p);
}

Vec CylindricalSurface::dfuv(const ParametricPoint& p) const
{
    return pimpl_->dfuv(p);
}

ostream& CylindricalSurface::print(ostream& os) const
{
    return pimpl_->print(os);
}

ParametricPoint CylindricalSurface::project(const Point& p) const
{
    return pimpl_->project(p);
}

CylindricalSurface::Impl::Impl()
    : r_(1)
    , ax_()
{
}

CylindricalSurface::Impl::Impl(double r, Axis ax)
    : r_(r)
    , ax_(move(ax))
{
}

Point CylindricalSurface::Impl::f(const ParametricPoint& p) const
{
    return r_ * ax_.pglobal(cos(p.u), sin(p.u), p.v);
}

Vec CylindricalSurface::Impl::dfu(const ParametricPoint& p) const
{
    return r_ * ax_.vglobal(-sin(p.u), cos(p.u), 0);
}

Vec CylindricalSurface::Impl::dfv(const ParametricPoint& p) const
{
    return ax_.vglobal(0, 0, r_);
}

Vec CylindricalSurface::Impl::dfuu(const ParametricPoint& p) const
{
    return -r_ * ax_.vglobal(cos(p.u), sin(p.u), 0);
}

Vec CylindricalSurface::Impl::dfuv(const ParametricPoint& p) const
{
    return ax_.vglobal(0, 0, 0);
}

Vec CylindricalSurface::Impl::dfvv(const ParametricPoint& p) const
{
    return ax_.vglobal(0, 0, 0);
}

ostream& CylindricalSurface::Impl::print(ostream& os) const
{
    fmt::print(os,
               "{{ \"type\": \"cylindrical\", \"r\": {0}, \"axis\": {1} }}",
               r_, ax_);
    return os;
}

ParametricPoint CylindricalSurface::Impl::project(const Point& p) const
{
    auto [c, x, y, z] = ax_.get_view();
    auto v = dot(p - c, z) / r_;
    auto w = (p - c) - z * dot(p - c, z);
    auto first = atan2(dot(w, y), dot(w, x));
    if (first < 0) {
        first += M_PI;
    }
    auto second = first + M_PI;

    return (dist(f({first, v}), p) < dist(f({second, v}), p))
        ? ParametricPoint{first, v}
        : ParametricPoint{second, v};
}