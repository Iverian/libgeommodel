#include <gm/edge.h>
#include <gm/face.h>
#include <gm/compare.h>

#include <util/debug.h>
#include <util/math.h>
#include <util/to_string.h>
#include <util/util.h>

#include <algorithm>
#include <stdexcept>

#include <fmt/ostream.h>

using namespace std;

namespace gm {

Edge::Edge(shared_ptr<AbstractCurve> curve, const Point& front,
           const Point& back)
    : pfront_()
    , pback_()
    , curve_(curve)
{
    pfront_ = curve_->project(front);
    pback_ = curve_->project_greater(back, pfront_).value();

    auto f = curve_->f(pfront_);
    auto b = curve_->f(pback_);

    if (!isnear(f, front, Tolerance::ZERO) || !isnear(b, back, Tolerance::ZERO)) {
        THROW_(runtime_error,
               "unable to construct edge with \"curve\": {0}, \"begin\": {1}, "
               "\"end\": {2}: projections \"begin_p\": {3}, \"end_p\": {4} do "
               "not match original points",
               *curve_, front, back, f, b);
    }
}

Edge::Edge(shared_ptr<AbstractCurve> curve, double pfront, double pback)
    : pfront_(pfront)
    , pback_(pback)
    , curve_(curve)
{
    if (pback_ < pfront_) {
        swap(pfront_, pback_);
    }
}

Point Edge::f(double u) const noexcept
{
    return curve_->f(param(u));
}

Vec Edge::df(double u) const noexcept
{
    return curve_->df(param(u));
}

Vec Edge::df2(double u) const noexcept
{
    return curve_->df2(param(u));
}

double Edge::project(const Point& p) const
{
    return param_rev(curve_->project(p));
}

optional<double> Edge::project_greater(const Point& p, double min) const
    noexcept
{
    auto result = curve_->project_greater(p, param(min));
    if (result.has_value()) {
        result = param_rev(result.value());
    }
    return result;
}

double Edge::approx_length(double begin, double end, size_t n) const
{
    return curve_->approx_length(param(begin), param(end), n);
}

Point Edge::front() const noexcept
{
    return curve_->f(pfront_);
}

Point Edge::back() const noexcept
{
    return curve_->f(pback_);
}

double Edge::pfront() const noexcept
{
    return pfront_;
}

double Edge::pback() const noexcept
{
    return pback_;
}

shared_ptr<AbstractCurve> Edge::curve() const noexcept
{
    return curve_;
}

bool operator==(const Edge& lhs, const Edge& rhs) noexcept
{
    return isnear(lhs.pfront_, rhs.pfront_) && isnear(lhs.pback_, rhs.pback_)
        && (lhs.curve_ == rhs.curve_);
}

bool operator!=(const Edge& lhs, const Edge& rhs) noexcept
{
    return !(lhs == rhs);
}

ostream& Edge::print(ostream& os) const
{
    fmt::print(os,
               "{{ \"curve\": {0}, \"front\": {1:.5g}, \"back\": {2:.5g} }}",
               *curve(), pfront_, pback_);
    return os;
}

double Edge::param(double t) const noexcept
{
    return pfront_ + t * (pback_ - pfront_);
}

double Edge::param_rev(double t) const noexcept
{
    return (t - pfront_) / (pback_ - pfront_);
}

} // namespace gm