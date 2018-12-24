#include <geom_model/line.h>
#include <util/debug.h>

#include <iostream>

#include <fmt/ostream.h>

using namespace std;

struct Line::Impl {
    Impl();
    Impl(Vec dir, Point c);

    const Point& center() const;
    const Vec& dir() const;

    Point f(double u) const;
    Vec df() const;
    Vec df2() const;
    double project(const Point& x) const;
    double length(double begin, double end) const;
    ostream& print(ostream& os) const;

private:
    Vec dir_;
    Point c_;
};

Line::~Line() = default;
Line::Line(Line&&) noexcept = default;
Line& Line::operator=(Line&&) noexcept = default;
Line::Line(const Line& other)
    : pimpl_(make_unique<Impl>(*other.pimpl_))
{
}
Line& Line::operator=(const Line& other)
{
    *pimpl_ = *other.pimpl_;
    return *this;
}

Line::Line()
    : pimpl_(make_unique<Impl>())
{
}

Line::Line(Vec dir, Point c)
    : pimpl_(make_unique<Impl>(move(dir), move(c)))
{
}

Point Line::f(double u) const
{
    return pimpl_->f(u);
}

Vec Line::df(double u) const
{
    return pimpl_->df();
}

Vec Line::df2(double u) const
{
    return pimpl_->df2();
}

double Line::project(const Point& p) const
{
    return pimpl_->project(p);
}

double Line::length(double begin, double end) const
{
    return pimpl_->length(begin, end);
}

ostream& Line::print(ostream& os) const
{
    return pimpl_->print(os);
}

const Point& Line::center() const
{
    return pimpl_->center();
}

const Vec& Line::dir() const
{
    return pimpl_->dir();
}

Line::Impl::Impl()
    : Impl({1, 0, 0}, {0, 0, 0})
{
}

Line::Impl::Impl(Vec dir, Point c)
    : dir_(move(dir))
    , c_(move(c))
{
}

const Point& Line::Impl::center() const
{
    return c_;
}

const Vec& Line::Impl::dir() const
{
    return dir_;
}

Point Line::Impl::f(double u) const
{
    return c_ + dir_ * u;
}

Vec Line::Impl::df() const
{
    return dir_;
}

Vec Line::Impl::df2() const
{
    return Vec();
}

ostream& Line::Impl::print(ostream& os) const
{

    fmt::print(os, "{{ \"type\": \"line\", \"dir\": {0}, \"center\": {1} }}",
               dir_, c_);
    return os;
}

double Line::Impl::project(const Point& p) const
{
    return dot(p - c_, Point(dir_)) / sqr(dir_);
}

double Line::Impl::length(double begin, double end) const
{
    return dir_.norm() * fabs(end - begin);
}

// vector<double> Line::Impl::discretize_param(double begin, double end) const
// {
//     static constexpr size_t mesh_size = 10;
//     vector<double> result;
//     for (size_t i = 0; i < mesh_size; ++i) {
//         result.emplace_back(begin
//                             + double(i) / (mesh_size - 1) * (end - begin));
//     }
//     return result;
// }