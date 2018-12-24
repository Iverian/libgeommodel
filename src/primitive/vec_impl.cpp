#include "vec_impl.h"
#include <util/debug.h>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <string>

EXCEPT(zero_division, "zero division")

using namespace std;

VecImpl::VecImpl()
    : data_{0, 0, 0}
{
}

VecImpl::VecImpl(double x, double y, double z)
    : data_{x, y, z}
{
}

VecImpl::VecImpl(initializer_list<double> data)
    : data_{}
{

    auto iter = data.begin();
    for (auto& i : data_)
        i = *iter++;
}

VecImpl::VecImpl(const VecImpl::data_type& data)
    : data_(data)
{
}

const double& VecImpl::operator[](size_t i) const noexcept
{
    return data_[i];
}
const double* VecImpl::data() const noexcept
{
    return data_.data();
}

size_t VecImpl::size() const noexcept
{
    return data_.size();
}

VecImpl& VecImpl::operator+=(const VecImpl& other)
{
    for (size_t i = 0; i < data_.size(); ++i)
        data_[i] += other.data_[i];
    return *this;
}

VecImpl& VecImpl::operator-=(const VecImpl& other)
{
    for (size_t i = 0; i < data_.size(); ++i)
        data_[i] -= other.data_[i];
    return *this;
}

VecImpl& VecImpl::operator*=(double x)
{
    for (auto& i : data_)
        i *= x;
    return *this;
}

VecImpl& VecImpl::operator/=(double x)
{
    return (*this *= 1 / x);
}

double VecImpl::dot(const VecImpl& other) const
{
    double result = 0;
    for (size_t i = 0; i < data_.size(); ++i)
        result += data_[i] * other.data_[i];
    return result;
}

VecImpl::data_type::const_iterator VecImpl::cbegin() const
{
    return data_.cbegin();
}

VecImpl::data_type::const_iterator VecImpl::cend() const
{
    return data_.cend();
}

VecImpl::data_type::iterator VecImpl::begin()
{
    return data_.begin();
}

VecImpl::data_type::iterator VecImpl::end()
{
    return data_.end();
}

bool operator==(const VecImpl& lhs, const VecImpl& rhs)
{
    return equal(lhs.cbegin(), lhs.cend(), rhs.cbegin());
}

bool operator!=(const VecImpl& lhs, const VecImpl& rhs)
{
    return !(lhs == rhs);
}

double VecImpl::sqr() const
{
    return dot(*this);
}

bool VecImpl::isnan() const
{
    return accumulate(
        cbegin(), cend(), false,
        [](const auto& lhs, const auto& rhs) { return lhs || ::isnan(rhs); });
}

const VecImpl::data_type& VecImpl::raw() const noexcept
{
    return data_;
}

double dot(const VecImpl& a, const VecImpl& b)
{
    return a.dot(b);
}

ostream& operator<<(ostream& os, const VecImpl& v)
{
    return os << "[" << v[0] << ", " << v[1] << ", " << v[2] << "]";
}

double sqr(const VecImpl& a)
{
    return a.sqr();
}
