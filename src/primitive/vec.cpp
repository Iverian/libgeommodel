#include <geom_model/point.h>
#include <geom_model/vec.h>

#include "vec_impl.h"
#include <util/math.h>

using namespace std;

Vec::Vec()
    : pimpl_(make_unique<VecImpl>())
{
}

Vec::Vec(double x, double y, double z)
    : pimpl_(make_unique<VecImpl>(x, y, z))
{
}

Vec::Vec(double magnitude, const array<double, 3>& coord)
    : Vec(coord)
{
    *this *= magnitude;
}

Vec::Vec(const initializer_list<double>& list)
    : pimpl_(make_unique<VecImpl>(list))
{
}

Vec::Vec(const array<double, 3>& coord)
    : pimpl_(make_unique<VecImpl>(coord))
{
}

Vec::Vec(const Point& a, const Point& b)
    : Vec(b - a)
{
}

Vec::Vec(const pair<Point, Point>& diff)
    : Vec(diff.second - diff.first)
{
}

Vec::Vec(const Point& p)
    : pimpl_(p.pimpl_)
{
}

Vec::Vec(Point&& p)
    : pimpl_(move(p.pimpl_))
{
}

array<double, 3>::iterator Vec::begin()
{
    copy();
    return pimpl_->begin();
}

array<double, 3>::iterator Vec::end()
{
    copy();
    return pimpl_->end();
}

array<double, 3>::const_iterator Vec::cbegin() const
{
    return pimpl_->cbegin();
}

array<double, 3>::const_iterator Vec::cend() const
{
    return pimpl_->cend();
}

Vec& Vec::operator+=(const Vec& other)
{
    copy();
    *pimpl_ += *other.pimpl_;
    return *this;
}

Vec& Vec::operator-=(const Vec& other)
{
    copy();
    *pimpl_ -= *other.pimpl_;
    return *this;
}

Vec& Vec::operator*=(double x)
{
    copy();
    *pimpl_ *= x;
    return *this;
}

Vec& Vec::operator/=(double x)
{
    copy();
    *pimpl_ /= x;
    return *this;
}

Vec Vec::operator+(const Vec& other) const
{
    Vec result = *this;
    return (result += other);
}

Vec Vec::operator-() const
{
    Vec result = *this;
    return (result *= -1);
}

Vec Vec::operator-(const Vec& other) const
{
    Vec result = *this;
    return (result -= other);
}

Vec Vec::operator*(double x) const
{
    Vec result = *this;
    return (result *= x);
}

Vec Vec::operator/(double x) const
{
    Vec result = *this;
    return (result /= x);
}

double Vec::sqr() const
{
    return pimpl_->sqr();
}

double Vec::norm() const
{
    return sqrt(sqr());
}

Vec& Vec::normalize()
{
    copy();
    return (*this /= norm());
}

Vec Vec::normalize() const
{
    Vec result = *this;
    return result.normalize();
}

bool operator==(const Vec& lhs, const Vec& rhs)
{
    return *lhs.pimpl_ == *rhs.pimpl_;
}

bool operator!=(const Vec& lhs, const Vec& rhs)
{
    return !(lhs == rhs);
}

const double& Vec::operator[](size_t i) const noexcept
{
    return (*pimpl_)[i];
}

size_t Vec::size() const noexcept
{
    return pimpl_->size();
}

const array<double, 3>& Vec::raw() const noexcept
{
    return pimpl_->raw();
}

ostream& operator<<(ostream& os, const Vec& v)
{
    return os << *v.pimpl_;
}

void Vec::copy()
{
    if (pimpl_.use_count() > 1) {
        pimpl_ = make_shared<VecImpl>(*pimpl_);
    }
}

double dist(const Vec& x, const Vec& y)
{
    return norm(x - y);
}

double norm(const Vec& x)
{
    return x.norm();
}

double sqr(const Vec& a)
{
    return a.sqr();
}

double cos(const Vec& a, const Vec& b)
{
    return dot(a, b) / (norm(a) * norm(b));
}

double sin(const Vec& a, const Vec& b)
{
    return norm(cross(a, b)) / (norm(a) * norm(b));
}

double angle(const Vec& a, const Vec& b)
{
    auto x = dot(a, b), y = norm(cross(a, b));
    return atan2(y, x);
}

bool Vec::isnan() const
{
    return pimpl_->isnan();
}

Vec cross(const Vec& a, const Vec& b)
{
    return Vec(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
               a[0] * b[1] - a[1] * b[0]);
}

Vec operator*(double x, const Vec& v)
{
    return v * x;
}

Vec unit(const Vec& v)
{
    return v.normalize();
}

Vec unit(Vec&& v)
{
    Vec result(move(v));
    result.normalize();
    return result;
}
