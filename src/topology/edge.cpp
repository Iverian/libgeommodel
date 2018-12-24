#include <geom_model/edge.h>
#include <geom_model/face.h>
#include <geom_model/geom_util.h>

#include <util/debug.h>
#include <util/math.h>
#include <util/to_string.h>
#include <util/util.h>

#include <algorithm>
#include <stdexcept>

#include <fmt/ostream.h>

using namespace std;

struct Edge::Impl {
    Impl();
    Impl(const shared_ptr<AbstractCurve>& curve, const Point& begin,
         const Point& end);
    Impl(const shared_ptr<AbstractCurve>& curve, double begin, double end);

    Point f(double t) const;
    Point front() const;
    Point back() const;
    double pfront() const;
    double pback() const;
    double length() const;

    const AbstractCurve& curve() const;

    bool operator==(const Impl& rhs) const;

private:
    double begin_;
    double end_;
    shared_ptr<AbstractCurve> curve_;
    mutable double len_;
};

Edge::Edge()
    : pimpl_(make_shared<Edge::Impl>())
{
}

Edge::Edge(const shared_ptr<AbstractCurve>& curve, const Point& begin,
           const Point& end)
    : pimpl_(make_shared<Edge::Impl>(curve, begin, end))
{
}

Edge::Edge(const shared_ptr<AbstractCurve>& curve, double begin, double end)
    : pimpl_(make_shared<Edge::Impl>(curve, begin, end))
{
}

Point Edge::f(double t) const
{
    return pimpl_->f(t);
}

Point Edge::front() const
{
    return pimpl_->front();
}

Point Edge::back() const
{
    return pimpl_->back();
}

double Edge::pfront() const
{
    return pimpl_->pfront();
}

double Edge::pback() const
{
    return pimpl_->pback();
}

double Edge::length() const
{
    return pimpl_->length();
}

const AbstractCurve& Edge::curve() const
{
    return pimpl_->curve();
}

bool operator==(const Edge& lhs, const Edge& rhs)
{
    return (*lhs.pimpl_) == (*rhs.pimpl_);
}

bool operator!=(const Edge& lhs, const Edge& rhs)
{
    return !(lhs == rhs);
}

ostream& operator<<(ostream& os, const Edge& edge)
{
    fmt::print(os, "{{ \"curve\": {0}, \"start\": {1}, \"end\": {2} }}",
               edge.curve(), edge.pfront(), edge.pback());
    return os;
}

Edge::Impl::Impl()
    : begin_(0)
    , end_(0)
    , curve_(nullptr)
{
}

Edge::Impl::Impl(const shared_ptr<AbstractCurve>& curve, const Point& begin,
                 const Point& end)
    : begin_(0)
    , end_(0)
    , curve_(curve)
{
    static constexpr auto eps = 1e-3;

    begin_ = curve_->project(begin);
    end_ = curve_->project_greater(end, begin_).value();

    auto b = curve_->f(begin_), e = curve_->f(end_);
    if (!isnear(b, begin, eps) || !isnear(e, end, eps)) {
        THROW_(runtime_error,
               "unable to construct edge with \"curve\": {0}, \"begin\": {1}, "
               "\"end\": {2}: projections \"begin_p\": {3}, \"end_p\": {4} do "
               "not match original points",
               *curve, begin, end, b, e);
    }
}

Edge::Impl::Impl(const shared_ptr<AbstractCurve>& curve, double begin,
                 double end)
    : begin_(begin)
    , end_(end)
    , len_(-1)
    , curve_(move(curve))
{
}

Point Edge::Impl::f(const double t) const
{
    return curve_->f(t);
}

Point Edge::Impl::front() const
{
    return f(begin_);
}

Point Edge::Impl::back() const
{
    return f(end_);
}

double Edge::Impl::pfront() const
{
    return begin_;
}

double Edge::Impl::pback() const
{
    return end_;
}

double Edge::Impl::length() const
{
    if (len_ < 0) {
        len_ = curve_->length(begin_, end_);
    }
    return len_;
}

const AbstractCurve& Edge::Impl::curve() const
{
    return *curve_;
}

bool Edge::Impl::operator==(const Edge::Impl& rhs) const
{
    return (begin_ == rhs.begin_) && (end_ == rhs.end_)
        && (curve_ == rhs.curve_);
}
