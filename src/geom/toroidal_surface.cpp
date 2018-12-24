#include <geom_model/circle.h>
#include <geom_model/toroidal_surface.h>

#include <util/math.h>

#include <fmt/ostream.h>

using namespace std;

struct ToroidalSurface::Impl {
    Impl();
    Impl(double r0, double r1, Axis ax);

    Point f(const ParametricPoint& p) const;
    Vec dfu(const ParametricPoint& p) const;
    Vec dfv(const ParametricPoint& p) const;
    Vec dfuu(const ParametricPoint& p) const;
    Vec dfuv(const ParametricPoint& p) const;
    Vec dfvv(const ParametricPoint& p) const;
    ostream& print(ostream& os) const;
    ParametricPoint project(const Point& p) const;

private:
    double r0_;
    double r1_;
    Axis ax_;
};

ToroidalSurface::~ToroidalSurface() = default;
ToroidalSurface::ToroidalSurface(ToroidalSurface&&) noexcept = default;
ToroidalSurface& ToroidalSurface::operator=(ToroidalSurface&&) noexcept
    = default;

ToroidalSurface::ToroidalSurface(const ToroidalSurface& other)
    : pimpl_(make_unique<Impl>(*other.pimpl_))
{
}

ToroidalSurface& ToroidalSurface::operator=(const ToroidalSurface& other)
{
    *pimpl_ = *other.pimpl_;
    return *this;
}

ToroidalSurface::ToroidalSurface()
    : pimpl_(make_unique<Impl>())
{
}

ToroidalSurface::ToroidalSurface(double r0, double r1, Axis ax)
    : pimpl_(make_unique<Impl>(r0, r1, move(ax)))
{
}

Point ToroidalSurface::f(const ParametricPoint& p) const
{
    return pimpl_->f(p);
}

Vec ToroidalSurface::dfu(const ParametricPoint& p) const
{
    return pimpl_->dfu(p);
}

Vec ToroidalSurface::dfv(const ParametricPoint& p) const
{
    return pimpl_->dfv(p);
}

Vec ToroidalSurface::dfuu(const ParametricPoint& p) const
{
    return pimpl_->dfuu(p);
}

Vec ToroidalSurface::dfvv(const ParametricPoint& p) const
{
    return pimpl_->dfvv(p);
}

Vec ToroidalSurface::dfuv(const ParametricPoint& p) const
{
    return pimpl_->dfuv(p);
}

ostream& ToroidalSurface::print(ostream& os) const
{
    return pimpl_->print(os);
}

ParametricPoint ToroidalSurface::project(const Point& p) const
{
    return pimpl_->project(p);
}

ToroidalSurface::Impl::Impl()
    : r0_(0.5)
    , r1_(1)
    , ax_()
{
}

ToroidalSurface::Impl::Impl(double r0, double r1, Axis ax)
    : r0_(r0)
    , r1_(r1)
    , ax_(move(ax))
{
}

Point ToroidalSurface::Impl::f(const ParametricPoint& p) const
{
    return ax_.pglobal((r1_ + r0_ * cos(p.v)) * cos(p.u),
                       (r1_ + r0_ * cos(p.v)) * sin(p.u), r0_ * sin(p.v));
}

Vec ToroidalSurface::Impl::dfu(const ParametricPoint& p) const
{
    return (r1_ + r0_ * cos(p.v)) * ax_.vglobal(-sin(p.u), cos(p.u), 0);
}

Vec ToroidalSurface::Impl::dfv(const ParametricPoint& p) const
{
    return r0_
        * ax_.vglobal(sin(p.v) * cos(p.u), sin(p.v) * sin(p.u), cos(p.v));
}

Vec ToroidalSurface::Impl::dfuu(const ParametricPoint& p) const
{
    return -(r1_ + r0_ * cos(p.v)) * ax_.vglobal(cos(p.u), sin(p.u), 0);
}

Vec ToroidalSurface::Impl::dfuv(const ParametricPoint& p) const
{
    return r0_ * sin(p.v) * ax_.vglobal(sin(p.u), cos(p.u), 0);
}

Vec ToroidalSurface::Impl::dfvv(const ParametricPoint& p) const
{
    return -r0_
        * ax_.vglobal(cos(p.v) * cos(p.u), cos(p.v) * sin(p.u), sin(p.v));
}

ostream& ToroidalSurface::Impl::print(ostream& os) const
{
    fmt::print(
        os,
        "{{ \"type\": \"toroidal\", \"r0\":{0}, \"r1\":{1}, \"axis\":{2} }}",
        r0_, r1_, ax_);
    return os;
}

ParametricPoint ToroidalSurface::Impl::project(const Point& p) const
{
    auto [c, x, y, z] = ax_.get_view();
    auto u = Circle((r0_ + r1_) / 2, ax_).project(p);
    auto v = Circle(r0_,
                    Axis::from_xy(cos(u) * x + sin(u) * y, z,
                                  c + r1_ * cos(u) * x + r1_ * sin(u) * y))
                 .project(p);
    return ParametricPoint{u, v};
}