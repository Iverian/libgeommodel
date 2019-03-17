#include <gm/compare.hpp>
#include <gm/edge.hpp>
#include <gm/face.hpp>

#include <util/debug.hpp>
#include <util/math.hpp>
#include <util/to_string.hpp>
#include <util/util.hpp>

#include <algorithm>
#include <stdexcept>

#include <fmt/ostream.h>

namespace gm {

Edge::Edge()
    : pfront_(0)
    , pback_(0)
    , curve_(0)
{
}

Edge::Edge(std::shared_ptr<AbstractCurve> curve, const Point& front,
           const Point& back)
    : pfront_()
    , pback_()
    , curve_(curve)
{
    pfront_ = curve_->project(front);
    pback_ = curve_->project_greater(back, pfront_).value();

    auto f = curve_->f(pfront_);
    auto b = curve_->f(pback_);

    check_if(
        cmp::near(f, front, Tolerance::ZERO)
            && cmp::near(b, back, Tolerance::ZERO),
        "unable to construct edge with \"curve\": {0}, \"std::begin\": {1}, "
        "\"std::end\": {2}: projections \"std::begin_p\": {3}, "
        "\"std::end_p\": {4} do "
        "not match original points",
        *curve_, front, back, f, b);
}

Edge::Edge(std::shared_ptr<AbstractCurve> curve, double pfront, double pback)
    : pfront_(pfront)
    , pback_(pback)
    , curve_(curve)
{
    if (pback_ < pfront_) {
        std::swap(pfront_, pback_);
    }
}

bool Edge::empty() const noexcept
{
    return curve_ == nullptr;
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

std::optional<double> Edge::project_greater(const Point& p, double min) const
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

std::shared_ptr<AbstractCurve> Edge::curve() const noexcept
{
    return curve_;
}

bool operator==(const Edge& lhs, const Edge& rhs) noexcept
{
    return cmp::near(lhs.pfront_, rhs.pfront_)
        && cmp::near(lhs.pback_, rhs.pback_) && (lhs.curve_ == rhs.curve_);
}

bool operator!=(const Edge& lhs, const Edge& rhs) noexcept
{
    return !(lhs == rhs);
}

std::ostream& Edge::print(std::ostream& os) const
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