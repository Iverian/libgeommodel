#include <geom_model/oriented_edge.h>

#include <fmt/ostream.h>

#include <util/to_string.h>

using namespace std;

struct OrientedEdge::Impl {
    Impl(Edge edge, bool orientation);
    const Edge& edge() const;
    bool orientation() const;
    Point f(double t) const;
    Point front() const;
    Point back() const;
    double pfront() const;
    double pback() const;
    double length() const;

private:
    Edge edge_;
    bool orientation_;

    mutable vector<Point> discrete_;
    mutable vector<double> discrete_param_;
};

OrientedEdge::OrientedEdge(const Edge& edge, bool orientation)
    : pimpl_(make_unique<Impl>(edge, orientation))
{
}

OrientedEdge::OrientedEdge(Edge&& edge, bool orienation)
    : pimpl_(make_unique<Impl>(move(edge), orienation))
{
}

const Edge& OrientedEdge::edge() const
{
    return pimpl_->edge();
}

bool OrientedEdge::orienation() const
{
    return pimpl_->orientation();
}

Point OrientedEdge::f(double t) const
{
    return pimpl_->f(t);
}

Point OrientedEdge::front() const
{
    return pimpl_->front();
}

Point OrientedEdge::back() const
{
    return pimpl_->back();
}

double OrientedEdge::pfront() const
{
    return pimpl_->pfront();
}

double OrientedEdge::pback() const
{
    return pimpl_->pback();
}

double OrientedEdge::length() const
{
    return pimpl_->length();
}

OrientedEdge::~OrientedEdge() = default;

OrientedEdge::OrientedEdge(OrientedEdge&&) noexcept = default;

OrientedEdge& OrientedEdge::operator=(OrientedEdge&&) noexcept = default;

OrientedEdge::OrientedEdge(const OrientedEdge& other)
    : pimpl_(make_unique<Impl>(*other.pimpl_))
{
}

OrientedEdge& OrientedEdge::operator=(const OrientedEdge& other)
{
    *pimpl_ = *other.pimpl_;
    return *this;
}

OrientedEdge::Impl::Impl(Edge edge, bool orientation)
    : edge_(move(edge))
    , orientation_(orientation)
{
}

const Edge& OrientedEdge::Impl::edge() const
{
    return edge_;
}

bool OrientedEdge::Impl::orientation() const
{
    return orientation_;
}

Point OrientedEdge::Impl::f(const double t) const
{
    return edge_.f(orientation_ ? t : 1 - t);
}

Point OrientedEdge::Impl::front() const
{
    return orientation_ ? edge_.front() : edge_.back();
}

Point OrientedEdge::Impl::back() const
{
    return orientation_ ? edge_.back() : edge_.front();
}

double OrientedEdge::Impl::pfront() const
{
    return orientation_ ? edge_.pfront() : edge_.pback();
}

double OrientedEdge::Impl::pback() const
{
    return orientation_ ? edge_.pback() : edge_.pfront();
}

double OrientedEdge::Impl::length() const
{
    return edge_.length();
}

ostream& operator<<(ostream& os, const OrientedEdge& oedge)
{
    fmt::print(os, "{{ \"edge\": {0}, \"orientation\": {1} }}", oedge.edge(),
               oedge.orienation());
    return os;
}