#include <geom_model/parabola.h>

#include <util/math.h>

#include <fmt/ostream.h>

using namespace std;

struct Parabola::Impl {
    Impl();
    Impl(double f, Axis ax);

    Point f(double u) const;
    Vec df(double u) const;
    Vec df2(double u) const;
    ostream& print(ostream& os) const;
    double project(const Point& p) const;
    double length(double begin, double end) const;

private:
    double f_;
    Axis ax_;
};

Parabola::~Parabola() = default;
Parabola::Parabola(Parabola&&) noexcept = default;
Parabola& Parabola::operator=(Parabola&&) noexcept = default;

Parabola::Parabola(const Parabola& other)
    : pimpl_(make_unique<Impl>(*other.pimpl_))
{
}

Parabola& Parabola::operator=(const Parabola& other)
{
    *pimpl_ = *other.pimpl_;
    return *this;
}

Parabola::Parabola()
    : pimpl_(make_unique<Impl>())
{
}

Parabola::Parabola(double f, Axis ax)
    : pimpl_(make_unique<Impl>(f, ax))
{
}

Point Parabola::f(double u) const
{
    return pimpl_->f(u);
}

Vec Parabola::df(double u) const
{
    return pimpl_->df(u);
}

Vec Parabola::df2(double u) const
{
    return pimpl_->df2(u);
}

ostream& Parabola::print(ostream& os) const
{
    return pimpl_->print(os);
}

double Parabola::project(const Point& p) const
{
    return pimpl_->project(p);
}

Parabola::Impl::Impl()
    : f_(1)
    , ax_()
{
}

Parabola::Impl::Impl(double f, Axis ax)
    : f_(fabs(f))
    , ax_(move(ax))
{
}

Point Parabola::Impl::f(double u) const
{
    return ax_.pglobal(f_ * sqr(u), 2 * u * f_, 0);
}

Vec Parabola::Impl::df(double u) const
{
    return 2 * f_ * ax_.vglobal(u, 1, 0);
}

Vec Parabola::Impl::df2(double u) const
{
    return ax_.vglobal(2 * f_, 0, 0);
}

ostream& Parabola::Impl::print(ostream& os) const
{
    fmt::print(os, "{{ \"type\": \"parabola\", \"f\": {0}, \"axis\": {1} }}", f_,
               ax_);
    return os;
}
// Проекция -- решение уравнения u^3 + a u - b == 0
double Parabola::Impl::project(const Point& p) const
{
    auto [c, x, y, z] = ax_.get_view();
    auto w = (p - c) - z * dot(p - c, z);
    auto a = (2 * f_ - dot(w, x)) / f_, b = dot(w, y) / f_;

    auto d = 9 * b + sqrt(12 * pow(a, 3) + 81 * sqr(b));
    auto numer = (-2 * pow(3, 1. / 3) * a + pow(2, 1. / 3) * pow(d, 2. / 3));
    auto denom = pow(6, 2. / 3) * pow(d, 1. / 3);
    return numer / denom;
}

#define sign(x) (signbit(x) ? (-1) : (1))

double Parabola::Impl::length(double begin, double end) const
{
    return abs(f_)
        * (pow(4 + sqr(begin), 1.5) / sign(begin) + pow(4 + sqr(end), 1.5),
           sign(end));
}

#undef sign