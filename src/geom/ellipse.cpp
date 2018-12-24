#include <geom_model/axis.h>
#include <geom_model/ellipse.h>

#include <util/math.h>

#include <fmt/ostream.h>

using namespace std;

struct Ellipse::Impl {
    Impl();
    Impl(double rx, double ry, Axis ax);

    double rx() const;
    double ry() const;
    const Axis& ax() const;

    Point f(double u) const;
    Vec df(double u) const;
    Vec df2(double u) const;
    double project(const Point& p) const;
    std::optional<double> project_greater(const Point& x, double min) const;
    ostream& print(ostream& os) const;

private:
    double rx_;
    double ry_;
    Axis ax_;
};

Ellipse::~Ellipse() = default;
Ellipse::Ellipse(Ellipse&&) noexcept = default;
Ellipse& Ellipse::operator=(Ellipse&&) noexcept = default;

Ellipse::Ellipse(const Ellipse& other)
    : pimpl_(make_unique<Impl>(*other.pimpl_)){};

Ellipse& Ellipse::operator=(const Ellipse& other)
{
    *pimpl_ = *other.pimpl_;
    return *this;
}

Ellipse::Ellipse()
    : pimpl_(make_unique<Impl>()){};
Ellipse::Ellipse(double rx, double ry, Axis ax)
    : pimpl_(make_unique<Impl>(rx, ry, ax)){};
;

double Ellipse::rx() const
{
    return pimpl_->rx();
}

double Ellipse::ry() const
{
    return pimpl_->ry();
}

const Axis& Ellipse::ax() const
{
    return pimpl_->ax();
}

Point Ellipse::f(double u) const
{
    return pimpl_->f(u);
};

Vec Ellipse::df(double u) const
{
    return pimpl_->df(u);
}

Vec Ellipse::df2(double u) const
{
    return pimpl_->df2(u);
}

double Ellipse::project(const Point& x) const
{
    return pimpl_->project(x);
}

optional<double> Ellipse::project_greater(const Point& x, double min) const
{
    return pimpl_->project_greater(x, min);
}

ostream& Ellipse::print(ostream& os) const
{
    return pimpl_->print(os);
}

Ellipse::Impl::Impl()
    : rx_(1)
    , ry_(1)
    , ax_()
{
}

Ellipse::Impl::Impl(double rx, double ry, Axis ax)
    : rx_(rx)
    , ry_(ry)
    , ax_(move(ax))
{
}

double Ellipse::Impl::rx() const
{
    return rx_;
}

double Ellipse::Impl::ry() const
{
    return ry_;
}

const Axis& Ellipse::Impl::ax() const
{
    return ax_;
}

Point Ellipse::Impl::f(double u) const
{
    return ax_.pglobal(rx_ * cos(u), ry_ * sin(u), 0);
}

Vec Ellipse::Impl::df(double u) const
{
    return ax_.vglobal(-rx_ * sin(u), ry_ * cos(u), 0);
}

Vec Ellipse::Impl::df2(double u) const
{
    return ax_.vglobal(-rx_ * cos(u), -ry_ * sin(u), 0);
}

ostream& Ellipse::Impl::print(ostream& os) const
{
    fmt::print(
        os,
        "{{ \"type\": \"ellipse\", \"rx\": {0}, \"ry\": {1}, \"axis\": {2} }}",
        rx_, ry_, ax_);
    return os;
}

double Ellipse::Impl::project(const Point& p) const
{
    auto x = Point(ax_[0]), y = Point(ax_[1]), z = Point(ax_[2]);
    auto c = ax_.center();
    auto w = (p - c) - z * dot(p - c, z);

    auto first = atan2(rx_ * dot(w, y), ry_ * dot(w, x));
    if (first < 0) {
        first += M_PI;
    }
    auto second = first + M_PI;

    return dist(f(first), p) < dist(f(second), p) ? first : second;
}

std::optional<double> Ellipse::Impl::project_greater(const Point& x,
                                                     double min) const
{
    auto result = project(x);
    while (result < min) {
        result += 2 * M_PI;
    }
    return result;
}