#include <geom_model/circle.h>

#include <util/math.h>

#include <fmt/ostream.h>

using namespace std;

struct Circle::Impl {
    double length(double begin, double end, const Ellipse& parent) const;
    ostream& print(ostream& os, const Ellipse& parent) const;
};

Circle::~Circle() = default;
Circle::Circle(Circle&&) noexcept = default;
Circle& Circle::operator=(Circle&&) noexcept = default;

Circle::Circle(const Circle& other)
    : Ellipse(other)
    , pimpl_(make_unique<Impl>(*other.pimpl_))
{
}

Circle& Circle::operator=(const Circle& other)
{
    Ellipse::operator=(other);
    *pimpl_ = *other.pimpl_;
    return *this;
}

Circle::Circle()
    : Ellipse()
    , pimpl_(make_unique<Impl>())
{
}

Circle::Circle(double r, Axis ax)
    : Ellipse(r, r, move(ax))
    , pimpl_(make_unique<Impl>())
{
}

double Circle::r() const
{
    return rx();
}

double Circle::length(double begin, double end) const
{
    return pimpl_->length(begin, end, *this);
}

ostream& Circle::print(ostream& os) const
{
    return pimpl_->print(os, *this);
}

double Circle::Impl::length(double begin, double end,
                            const Ellipse& parent) const
{
    return parent.rx() * fabs(end - begin);
}

ostream& Circle::Impl::print(ostream& os, const Ellipse& parent) const
{
    fmt::print(os, "{{ \"type\": \"circle\", \"r\": {0}, \"axis\": {1} }}",
               parent.rx(), parent.ax());
    return os;
}
