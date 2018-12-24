#include <geom_model/plane.h>

#include <fmt/ostream.h>

using namespace std;

struct Plane::Impl {
    Impl();
    explicit Impl(Axis ax);

    const Axis& ax() const;

    Point f(const ParametricPoint& p) const;
    Vec dfu(const ParametricPoint& p) const;
    Vec dfv(const ParametricPoint& p) const;
    Vec dfuu(const ParametricPoint& p) const;
    Vec dfuv(const ParametricPoint& p) const;
    Vec dfvv(const ParametricPoint& p) const;
    ostream& print(ostream& os) const;
    ParametricPoint project(const Point& p) const;

private:
    Axis ax_;
};

Plane::~Plane() = default;
Plane::Plane(Plane&&) noexcept = default;
Plane& Plane::operator=(Plane&&) noexcept = default;

Plane::Plane(const Plane& other)
    : pimpl_(make_unique<Impl>(*other.pimpl_))
{
}

Plane& Plane::operator=(const Plane& other)
{
    *pimpl_ = *other.pimpl_;
    return *this;
}

Plane::Plane()
    : pimpl_(make_unique<Impl>())
{
}

Plane::Plane(Axis ax)
    : pimpl_(make_unique<Impl>(move(ax)))
{
}

const Axis& Plane::ax() const
{
    return pimpl_->ax();
}

Point Plane::f(const ParametricPoint& p) const
{
    return pimpl_->f(p);
}

Vec Plane::dfu(const ParametricPoint& p) const
{
    return pimpl_->dfu(p);
}

Vec Plane::dfv(const ParametricPoint& p) const
{
    return pimpl_->dfv(p);
}

Vec Plane::dfuu(const ParametricPoint& p) const
{
    return pimpl_->dfuu(p);
}

Vec Plane::dfvv(const ParametricPoint& p) const
{
    return pimpl_->dfvv(p);
}

Vec Plane::dfuv(const ParametricPoint& p) const
{
    return pimpl_->dfuv(p);
}

ostream& Plane::print(ostream& os) const
{
    return pimpl_->print(os);
}

ParametricPoint Plane::project(const Point& p) const
{
    return pimpl_->project(p);
}

Plane::Impl::Impl()
    : ax_()
{
}

Plane::Impl::Impl(Axis ax)
    : ax_(move(ax))
{
}

const Axis& Plane::Impl::ax() const
{
    return ax_;
}

Point Plane::Impl::f(const ParametricPoint& p) const
{
    return ax_.pglobal(p.u, p.v, 0);
}

Vec Plane::Impl::dfu(const ParametricPoint& p) const
{
    return ax_.vglobal(1, 0, 0);
}

Vec Plane::Impl::dfv(const ParametricPoint& p) const
{
    return ax_.vglobal(0, 1, 0);
}

Vec Plane::Impl::dfuu(const ParametricPoint& p) const
{
    return ax_.vglobal(0, 0, 0);
}

Vec Plane::Impl::dfuv(const ParametricPoint& p) const
{
    return ax_.vglobal(0, 0, 0);
}

Vec Plane::Impl::dfvv(const ParametricPoint& p) const
{
    return ax_.vglobal(0, 0, 0);
}

ostream& Plane::Impl::print(ostream& os) const
{
    fmt::print(os, "{{ \"type\": \"plane\", \"axis\": {0} }}", ax_);
    return os;
}

ParametricPoint Plane::Impl::project(const Point& p) const
{
    auto [c, x, y, z] = ax_.get_view();
    auto w = (p - c) - z * dot(p - c, z);
    return {dot(w, x), dot(w, y)};
}