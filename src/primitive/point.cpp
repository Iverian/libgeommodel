#include <geom_model/point.h>
#include <geom_model/vec.h>

#include "vec_impl.h"
#include <util/math.h>

#include <algorithm>
#include <array>

using namespace std;

Point::Point()
    : pimpl_(make_unique<VecImpl>())
{
}

Point::Point(double x, double y, double z)
    : pimpl_(make_unique<VecImpl>(x, y, z))
{
}

Point::Point(const initializer_list<double>& list)
    : pimpl_(make_unique<VecImpl>(list))
{
}

Point::Point(const array<double, 3>& coord)
    : pimpl_(make_unique<VecImpl>(coord))
{
}

Point::Point(const Vec& v)
    : pimpl_(v.pimpl_)
{
}

Point::Point(Vec&& v)
    : pimpl_(move(v.pimpl_))
{
}

array<double, 3>::iterator Point::begin()
{
    copy();
    return pimpl_->begin();
}

array<double, 3>::iterator Point::end()
{
    copy();
    return pimpl_->end();
}

array<double, 3>::const_iterator Point::cbegin() const
{
    return pimpl_->cbegin();
}

array<double, 3>::const_iterator Point::cend() const
{
    return pimpl_->cend();
}

Point& Point::operator+=(const Point& other)
{
    copy();
    *pimpl_ += *other.pimpl_;
    return *this;
}

Point& Point::operator-=(const Point& other)
{
    copy();
    *pimpl_ -= *other.pimpl_;
    return *this;
}

Point& Point::operator*=(double x)
{
    copy();
    *pimpl_ *= x;
    return *this;
}

Point& Point::operator/=(double x)
{
    copy();
    *pimpl_ /= x;
    return *this;
}

Point Point::operator+(const Point& other) const
{
    Point result = *this;
    return (result += other);
}

Point Point::operator-() const
{
    Point result = *this;
    return (result *= -1);
}

Point Point::operator-(const Point& other) const
{
    Point result = *this;
    return (result -= other);
}

Point Point::operator*(double x) const
{
    Point result = *this;
    return (result *= x);
}

Point Point::operator/(double x) const
{
    Point result = *this;
    return (result /= x);
}

double Point::dist(const Point& other) const
{
    return sqrt((*this - other).pimpl_->sqr());
}

bool Point::isnear(const Point& other, double eps) const
{
    return dist(other) < eps;
}

bool Point::isnan() const
{
    return pimpl_->isnan();
}

bool operator==(const Point& lhs, const Point& rhs)
{
    return *lhs.pimpl_ == *rhs.pimpl_;
}

bool operator!=(const Point& lhs, const Point& rhs)
{
    return !(lhs == rhs);
}

const double& Point::operator[](size_t i) const noexcept
{
    return (*pimpl_)[i];
}

size_t Point::size() const noexcept
{
    return pimpl_->size();
}

const array<double, 3>& Point::raw() const noexcept
{
    return pimpl_->raw();
}

double dist(const Point& a, const Point& b)
{
    return a.dist(b);
}

double angle(const Point& a, const Point& b, const Point& c, const Vec& norm)
{
    auto ba = Vec(b, a), bc = Vec(b, c);
    auto result = angle(ba, bc);
    if (dot(cross(ba, bc), norm) > 0) {
        result = 2 * M_PI - result;
    }
    return result;
}

Point operator+(const Point& p, const Vec& v)
{
    return p + Point(v);
}

Point operator+(const Vec& v, const Point& p)
{
    return p + v;
}

Point operator-(const Point& p, const Vec& v)
{
    return p - Point(v);
}

Point operator-(const Vec& v, const Point& p)
{
    return Point(v) - p;
}

Point operator*(double x, const Point& p)
{
    return p * x;
}

ostream& operator<<(ostream& os, const Point& p)
{
    return os << *p.pimpl_;
}

bool isnear(const Point& lhs, const Point& rhs, double eps)
{
    return lhs.isnear(rhs, eps);
}

void Point::copy()
{
    if (pimpl_.use_count() > 1)
        pimpl_ = make_shared<VecImpl>(*pimpl_);
}
