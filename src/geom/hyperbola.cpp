#include <geom_model/hyperbola.h>

#include <util/math.h>

#include <fmt/ostream.h>

using namespace std;

struct Hyperbola::Impl {
    Impl();
    Impl(double rx, double ry, Axis ax);

    Point f(double u) const;
    Vec df(double u) const;
    Vec df2(double u) const;
    ostream& print(ostream& os) const;
    double project(const Point& p) const;

private:
    double rx_;
    double ry_;
    Axis ax_;
};

Hyperbola::~Hyperbola() = default;
Hyperbola::Hyperbola(Hyperbola&&) noexcept = default;
Hyperbola& Hyperbola::operator=(Hyperbola&&) noexcept = default;

Hyperbola::Hyperbola(const Hyperbola& other)
    : pimpl_(make_unique<Impl>(*other.pimpl_))
{
}

Hyperbola& Hyperbola::operator=(const Hyperbola& other)
{
    *pimpl_ = *other.pimpl_;
    return *this;
}

Hyperbola::Hyperbola()
    : pimpl_(make_unique<Impl>())
{
}

Hyperbola::Hyperbola(double rx, double ry, Axis ax)
    : pimpl_(make_unique<Impl>(rx, ry, move(ax)))
{
}

Point Hyperbola::f(double u) const
{
    return pimpl_->f(u);
}

Vec Hyperbola::df(double u) const
{
    return pimpl_->df(u);
}

Vec Hyperbola::df2(double u) const
{
    return pimpl_->df2(u);
}

ostream& Hyperbola::print(ostream& os) const
{
    return pimpl_->print(os);
}

double Hyperbola::project(const Point& p) const
{
    return pimpl_->project(p);
}

Hyperbola::Impl::Impl()
    : rx_(1)
    , ry_(1)
    , ax_()
{
}

Hyperbola::Impl::Impl(double rx, double ry, Axis ax)
    : rx_(rx)
    , ry_(ry)
    , ax_(move(ax))
{
}

Point Hyperbola::Impl::f(double u) const
{
    return ax_.pglobal(rx_ * cosh(u), ry_ * sinh(u), 0);
}

Vec Hyperbola::Impl::df(double u) const
{
    return ax_.vglobal(rx_ * sinh(u), ry_ * cosh(u), 0);
}

Vec Hyperbola::Impl::df2(double u) const
{
    return ax_.vglobal(rx_ * cosh(u), ry_ * sinh(u), 0);
}

ostream& Hyperbola::Impl::print(ostream& os) const
{
    fmt::print(os,
               "{{ \"type\": \"hyperbola\", \"rx\": {0}, \"ry\": {1}, "
               "\"axis\": {2} }}",
               rx_, ry_, ax_);
    return os;
}

double Hyperbola::Impl::project(const Point& p) const
{
    auto [c, x, y, z] = ax_.get_view();
    auto w = (p - c) - z * dot(p - c, z);
    return atanh((rx_ * dot(w, y)) / (ry_ * dot(w, x)));
}
