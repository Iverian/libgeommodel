#ifndef GEOM_MODEL_BSPLINE_WPOINT_H_
#define GEOM_MODEL_BSPLINE_WPOINT_H_

#include <bspline/array_operators.h>
#include <gm/point.h>
#include <util/debug.h>
#include <util/itertools.h>

#include <array>
#include <ostream>
#include <vector>

namespace gm {

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

    using Proj = std::array<T, N>;

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

    WPoint(const Proj& p, value_type w)
        : Super()
    {
        for (size_type i = 0; i < N; ++i) {
            (*this)[i] = p[i] * w;
        }
        (*this)[N] = w;
    }

    WPoint(const WPoint&) = default;
    WPoint(WPoint&&) noexcept = default;
    WPoint& operator=(const WPoint&) = default;
    WPoint& operator=(WPoint&&) noexcept = default;

    Proj wp() const noexcept
    {
        Proj result;
        for (size_type i = 0; i < N; ++i) {
            result[i] = (*this)[i];
        }
        return result;
    }

    Proj p() const __GM_NOEXCEPT_RELEASE__
    {
        auto& ww = w();
        check_ifd(!gm::cmp::zero(ww), "Zero weight");

        Proj result;
        for (size_type i = 0; i < N; ++i) {
            result[i] = (*this)[i] / ww;
        }
        return result;
    }

    const_reference w() const noexcept
    {
        return (*this)[N];
    }

    constexpr size_type size() const noexcept
    {
        return N;
    }

    WPoint& operator+=(const WPoint& rhs) noexcept
    {
        for (size_type i = 0; i < N + 1; ++i) {
            (*this)[i] += rhs[i];
        }
        return *this;
    }

    WPoint& operator-=(const WPoint& rhs) noexcept
    {
        for (size_type i = 0; i < N + 1; ++i) {
            (*this)[i] -= rhs[i];
        }
        return *this;
    }

    WPoint& operator*=(value_type rhs) noexcept
    {
        for (size_type i = 0; i < N + 1; ++i) {
            (*this)[i] *= rhs;
        }
        return *this;
    }

    WPoint& operator/=(value_type rhs) __GM_NOEXCEPT_RELEASE__
    {
        check_ifd(!gm::cmp::zero(rhs), "Division by zero");

        for (size_type i = 0; i < N + 1; ++i) {
            (*this)[i] /= rhs;
        }
        return *this;
    }
};

using CPoint = WPoint<double, 3>;
using DPoint = WPoint<double, 1>;

Point pget(const CPoint& p) __GM_NOEXCEPT_RELEASE__;
double pget(const DPoint& p) __GM_NOEXCEPT_RELEASE__;

template <class T, size_t N>
DPoint wdot(const WPoint<T, N>& lhs, const WPoint<T, N>& rhs) noexcept
{
    T wp {};
    for (size_t i = 0; i < N; ++i) {
        wp += lhs[i] * rhs[i];
    }
    return DPoint({wp, lhs[N] * rhs[N]});
}

template <class T, size_t N>
WPoint<T, N> operator+(const WPoint<T, N>& lhs,
                       const typename WPoint<T, N>::Proj& rhs)
{
    auto result = lhs;
    auto w = lhs.w();
    for (size_t i = 0; i < N; ++i) {
        result[i] += w * rhs[i];
    }
    return result;
}

template <class T, size_t N>
WPoint<T, N> operator+(const typename WPoint<T, N>::Proj& lhs,
                       const WPoint<T, N>& rhs)
{
    return rhs + lhs;
}

template <class T, size_t N>
WPoint<T, N> operator-(const WPoint<T, N>& lhs,
                       const typename WPoint<T, N>::Proj& rhs)
{
    auto result = lhs;
    auto w = lhs.w();
    for (size_t i = 0; i < N; ++i) {
        result[i] -= w * rhs[i];
    }
    return result;
}

template <class T, size_t N>
WPoint<T, N> operator-(const typename WPoint<T, N>::Proj& lhs,
                       const WPoint<T, N>& rhs)
{
    WPoint<T, N> result;
    auto w = rhs.w();
    for (size_t i = 0; i < N; ++i) {
        result[i] = lhs[i] * w - rhs[i];
    }
    result[N] = w;

    return result;
}

template <class T, size_t N>
WPoint<T, N> operator+(const WPoint<T, N>& lhs,
                       const WPoint<T, N>& rhs) noexcept
{
    auto result = lhs;
    return (result += rhs);
}

template <class T, size_t N>
WPoint<T, N> operator-(const WPoint<T, N>& lhs,
                       const WPoint<T, N>& rhs) noexcept
{
    auto result = lhs;
    return (result -= rhs);
}

template <class T, size_t N>
WPoint<T, N> operator*(const WPoint<T, N>& lhs,
                       typename WPoint<T, N>::const_reference rhs) noexcept
{
    auto result = lhs;
    return (result *= rhs);
}

template <class T, size_t N>
WPoint<T, N>
operator/(const WPoint<T, N>& lhs,
          typename WPoint<T, N>::const_reference rhs) __GM_NOEXCEPT_RELEASE__
{
    auto result = lhs;
    return (result /= rhs);
}

template <class T, size_t N>
WPoint<T, N> operator*(typename WPoint<T, N>::const_reference lhs,
                       const WPoint<T, N>& rhs) noexcept
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
[[nodiscard]] WPoint<T, N> diff1(const std::vector<WPoint<T, N>>& p) noexcept
{
    auto pp = p[1].wp() / p[0].w() - p[1].w() * p[0].wp() / sqr(p[0].w());
    return WPoint<T, N>(pp, 1);
}

template <class T, size_t N>
[[nodiscard]] WPoint<T, N> diff2(const std::vector<WPoint<T, N>>& p) noexcept
{
    auto pp = (p[1].wp() * (-2 * p[1].w() / p[0].w())
               + p[0].wp() * (2 * p[1].w() / sqr(p[0].w())) + p[2].wp()
               - p[0].wp() * (p[2].w() / p[0].w()))
        / p[0].w();
    return WPoint<T, N>(pp, 1);
}

} // namespace gm

#endif // GEOM_MODEL_BSPLINE_WPOINT_H_