#include <geom_model/conical_surface.h>
#include <geom_model/line.h>

#include <util/math.h>

#include <fmt/ostream.h>

using namespace std;

struct ConicalSurface::Impl {
    Impl();
    Impl(double r, double a, Axis ax);

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
    double ta_;
    Axis ax_;
};

ConicalSurface::~ConicalSurface() = default;
ConicalSurface::ConicalSurface(ConicalSurface&&) noexcept = default;
ConicalSurface& ConicalSurface::operator=(ConicalSurface&&) noexcept = default;

ConicalSurface::ConicalSurface(const ConicalSurface& other)
    : pimpl_(make_unique<Impl>(*other.pimpl_))
{
}

ConicalSurface& ConicalSurface::operator=(const ConicalSurface& other)
{
    *pimpl_ = *other.pimpl_;
    return *this;
}

ConicalSurface::ConicalSurface()
    : pimpl_(make_unique<Impl>())
{
}

ConicalSurface::ConicalSurface(double r, double a, Axis ax)
    : pimpl_(make_unique<Impl>(r, a, move(ax)))
{
}

Point ConicalSurface::f(const ParametricPoint& p) const
{
    return pimpl_->f(p);
}

Vec ConicalSurface::dfu(const ParametricPoint& p) const
{
    return pimpl_->dfu(p);
}

Vec ConicalSurface::dfv(const ParametricPoint& p) const
{
    return pimpl_->dfv(p);
}

Vec ConicalSurface::dfuu(const ParametricPoint& p) const
{
    return pimpl_->dfuu(p);
}

Vec ConicalSurface::dfvv(const ParametricPoint& p) const
{
    return pimpl_->dfvv(p);
}

Vec ConicalSurface::dfuv(const ParametricPoint& p) const
{
    return pimpl_->dfuv(p);
}

ostream& ConicalSurface::print(ostream& os) const
{
    return pimpl_->print(os);
}

ParametricPoint ConicalSurface::project(const Point& p) const
{
    return pimpl_->project(p);
}

ConicalSurface::Impl::Impl()
    : r_(1)
    , ta_(tan(M_PI / 6))
    , ax_()
{
}

ConicalSurface::Impl::Impl(double r, double a, Axis ax)
    : r_(r)
    , ta_(tan(a))
    , ax_(move(ax))
{
}

Point ConicalSurface::Impl::f(const ParametricPoint& p) const
{
    return ax_.pglobal((r_ + p.v * ta_) * cos(p.u),
                       (r_ + p.v * ta_) * sin(p.u), p.v);
}

Vec ConicalSurface::Impl::dfu(const ParametricPoint& p) const
{
    return ax_.vglobal(-(r_ + p.v * ta_) * sin(p.u),
                       (r_ + p.v * ta_) * cos(p.u), 0);
}

Vec ConicalSurface::Impl::dfv(const ParametricPoint& p) const
{
    return ax_.vglobal(ta_ * cos(p.u), ta_ * sin(p.u), 1);
}

Vec ConicalSurface::Impl::dfuu(const ParametricPoint& p) const
{
    return ax_.vglobal(-(r_ + p.v * ta_) * cos(p.u),
                       (r_ + p.v * ta_) * sin(p.u), 0);
}

Vec ConicalSurface::Impl::dfuv(const ParametricPoint& p) const
{
    return ax_.vglobal(-ta_ * sin(p.u), ta_ * cos(p.u), 0);
}

Vec ConicalSurface::Impl::dfvv(const ParametricPoint& p) const
{
    return ax_.vglobal(0, 0, 0);
}

ostream& ConicalSurface::Impl::print(ostream& os) const
{
    fmt::print(
        os,
        "{{ \"type\": \"conical\", \"r\": {0}, \"ta\": {1}, \"axis\": {2} }}",
        r_, ta_, ax_);
    return os;
}

ParametricPoint ConicalSurface::Impl::project(const Point& p) const
{
    auto [c, x, y, z] = ax_.get_view();
    auto w = (p - c) - z * dot(p - c, z);
    auto u = atan2v(dot(w, y), dot(w, x));
    array<double, 2> v;
    transform(cbegin(u), cend(u), begin(v), [&](const auto& i) {
        return Line(ta_ * cos(i) * x + ta_ * sin(i) * y + z,
                    r_ * cos(i) * x + r_ * sin(i) * y + c)
            .project(p);
    });

    return dist(f({u[0], v[0]}), p) < dist(f({u[1], v[1]}), p)
        ? ParametricPoint{u[0], v[0]}
        : ParametricPoint{u[1], v[1]};
}