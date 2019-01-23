#ifndef GEOM_MODEL_SRC_PRIMITIVE_WPOINT_H_
#define GEOM_MODEL_SRC_PRIMITIVE_WPOINT_H_

#include <array>
#include <ostream>
#include <vector>

#include <util/itertools.h>

template <class T, size_t N>
class WPoint : std::array<T, N + 1> {
    using Super = std::array<T, N + 1>;

public:
    using typename Super::const_pointer;
    using typename Super::const_reference;
    using typename Super::difference_type;
    using typename Super::pointer;
    using typename Super::reference;
    using typename Super::size_type;
    using typename Super::value_type;

    using Super::operator[];
    using Super::begin;
    using Super::end;
    using Super::size;

    WPoint()
        : Super()
    {
    }
    WPoint(const Super& value)
        : Super(value)
    {
    }
    WPoint(Super&& value)
        : Super(value)
    {
    }
    WPoint(const std::array<T, N>& p, value_type w)
        : Super()
    {
        for (size_type i = 0; i < N; ++i) {
            (*this)[i] = p[i] * w;
        }
        (*this)[N] = w;
    }

    std::array<T, N> wp() const
    {
        std::array<T, N> result;
        for (size_type i = 0; i < N; ++i) {
            result[i] = (*this)[i];
        }
        return result;
    }
    std::array<T, N> p() const
    {
        std::array<T, N> result;
        for (size_type i = 0; i < N; ++i) {
            result[i] = (*this)[i] / (*this)[N];
        }
        return result;
    }
    const_reference w() const
    {
        return (*this)[N];
    }

    WPoint& operator+=(const WPoint& rhs)
    {
        for (size_type i = 0; i < N + 1; ++i) {
            (*this)[i] += rhs[i];
        }
        return *this;
    }
    WPoint& operator-=(const WPoint& rhs)
    {
        for (size_type i = 0; i < N + 1; ++i) {
            (*this)[i] -= rhs[i];
        }
        return *this;
    }
    WPoint& operator*=(value_type rhs)
    {
        for (size_type i = 0; i < N + 1; ++i) {
            (*this)[i] *= rhs;
        }
        return *this;
    }
    WPoint& operator/=(value_type rhs)
    {
        for (size_type i = 0; i < N + 1; ++i) {
            (*this)[i] /= rhs;
        }
        return *this;
    }

    static WPoint d1(const std::vector<WPoint>& p)
    {
        auto result
            = p[1].wp() / p[0].w() - p[1].w() * p[0].wp() / sqr(p[0].w());
        return WPoint(result, 1);
    }
    static WPoint d2(const std::vector<WPoint>& p)
    {
        auto result = (p[1].wp() * (-2 * p[1].w() / p[0].w())
                       + p[0].wp() * (2 * p[1].w() / sqr(p[0].w())) + p[2].wp()
                       - p[0].wp() * (p[2].w() / p[0].w()))
            / p[0].w();
        return WPoint(result, 1);
    }
};

using CPoint = WPoint<double, 3>;
using DPoint = WPoint<double, 1>;

template <class T, size_t N>
WPoint<T, N> operator+(const WPoint<T, N>& lhs, const WPoint<T, N>& rhs)
{
    auto result = lhs;
    return (result += rhs);
}

template <class T, size_t N>
WPoint<T, N> operator-(const WPoint<T, N>& lhs, const WPoint<T, N>& rhs)
{
    auto result = lhs;
    return (result -= rhs);
}

template <class T, size_t N>
WPoint<T, N> operator*(const WPoint<T, N>& lhs,
                       typename WPoint<T, N>::const_reference rhs)
{
    auto result = lhs;
    return (result *= rhs);
}

template <class T, size_t N>
WPoint<T, N> operator/(const WPoint<T, N>& lhs,
                       typename WPoint<T, N>::const_reference rhs)
{
    auto result = lhs;
    return (result /= rhs);
}

template <class T, size_t N>
WPoint<T, N> operator*(typename WPoint<T, N>::const_reference lhs,
                       const WPoint<T, N>& rhs)
{
    auto result = rhs;
    return (result *= lhs);
}

template <class T, size_t N>
std::ostream& operator<<(std::ostream& os, const WPoint<T, N>& obj)
{
    fmt::print(os, "{{ \"wp\": {}, \"w\": {} }}",
               RangePrint(obj.begin(), std::prev(obj.end())), obj[N]);
    return os;
}

template <class T, size_t N>
std::array<T, N>& operator+=(std::array<T, N>& lhs,
                             const std::array<T, N>& rhs)
{
    for (typename std::array<T, N>::size_type i = 0; i < N; ++i) {
        lhs[i] += rhs[i];
    }
    return lhs;
}

template <class T, size_t N>
std::array<T, N>& operator-=(std::array<T, N>& lhs,
                             const std::array<T, N>& rhs)
{
    for (typename std::array<T, N>::size_type i = 0; i < N; ++i) {
        lhs[i] -= rhs[i];
    }
    return lhs;
}

template <class T, size_t N>
std::array<T, N>& operator*=(std::array<T, N>& lhs, const T& rhs)
{
    for (typename std::array<T, N>::size_type i = 0; i < N; ++i) {
        lhs[i] *= rhs;
    }
    return lhs;
}

template <class T, size_t N>
std::array<T, N>& operator/=(std::array<T, N>& lhs, const T& rhs)
{
    for (typename std::array<T, N>::size_type i = 0; i < N; ++i) {
        lhs[i] /= rhs;
    }
    return lhs;
}

template <class T, size_t N>
std::array<T, N> operator+(std::array<T, N> lhs, std::array<T, N> rhs)
{
    auto result = std::move(lhs);
    return (result += rhs);
}

template <class T, size_t N>
std::array<T, N> operator-(std::array<T, N> lhs, std::array<T, N> rhs)
{
    auto result = std::move(lhs);
    return (result -= rhs);
}

template <class T, size_t N>
std::array<T, N> operator*(std::array<T, N> lhs, T rhs)
{
    auto result = std::move(lhs);
    return (result *= rhs);
}

template <class T, size_t N>
std::array<T, N> operator*(T lhs, std::array<T, N> rhs)
{
    auto result = std::move(rhs);
    return (result *= lhs);
}

template <class T, size_t N>
std::array<T, N> operator/(std::array<T, N> lhs, T rhs)
{
    auto result = std::move(lhs);
    return (result /= rhs);
}

#endif // GEOM_MODEL_SRC_PRIMITIVE_WPOINT_H_