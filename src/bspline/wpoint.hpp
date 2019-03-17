#ifndef GEOM_MODEL_SRC_BSPLINE_WPOINT_HPP_
#define GEOM_MODEL_SRC_BSPLINE_WPOINT_HPP_

// #include <bspline/array_operators.h>
#include <util/debug.h>
#include <util/itertools.h>

#include <array>
#include <ostream>
#include <type_traits>
#include <vector>

namespace gm {

template <class T, size_t N>
class WPoint {
    using Data = std::array<T, N + 1>;

    Data d_;

public:
    using size_type = typename Data::size_type;
    using value_type = typename Data::value_type;
    using reference = typename Data::reference;
    using const_reference = typename Data::const_reference;
    using pointer = typename Data::pointer;
    using const_pointer = typename Data::const_pointer;
    using iterator = typename Data::iterator;
    using const_iterator = typename Data::const_iterator;
    // using typename Super::const_pointer;
    // using typename Super::const_reference;
    // using typename Super::difference_type;
    // using typename Super::pointer;
    // using typename Super::reference;
    // using typename Super::size_type;
    // using typename Super::value_type;

    // using Super::operator[];
    // using Super::begin;
    // using Super::end;
    // using Super::size;

    using Proj = std::array<T, N>;

    WPoint()
        : d_()
    {
    }

    WPoint(const Data& value)
        : d_(value)
    {
    }

    WPoint(Data&& value)
        : d_(value)
    {
    }

    WPoint(const Proj& p, value_type w)
        : d_()
    {
        for (size_type i = 0; i < N; ++i) {
            d_[i] = p[i] * w;
        }
        d_[N] = w;
    }

    WPoint(const WPoint&) = default;
    WPoint(WPoint&&) noexcept = default;
    WPoint& operator=(const WPoint&) = default;
    WPoint& operator=(WPoint&&) noexcept = default;

    reference operator[](size_type i) noexcept
    {
        return d_[i];
    }

    const_reference operator[](size_type i) const noexcept
    {
        return d_[i];
    }

    iterator begin() noexcept
    {
        return d_.begin();
    }

    iterator end() noexcept
    {
        return d_.end();
    }

    const_iterator begin() const noexcept
    {
        return d_.begin();
    }

    const_iterator end() const noexcept
    {
        return d_.end();
    }

    constexpr size_type size() const noexcept
    {
        return N + 1;
    }

    Proj wp() const noexcept
    {
        Proj result;
        for (size_type i = 0; i < N; ++i) {
            result[i] = d_[i];
        }
        return result;
    }

    WPoint<T, N>& wdiv(const_reference v) __GM_NOEXCEPT_RELEASE__
    {
        check_ifd(!gm::cmp::zero(v), "Division by zero");
        for (size_type i = 0; i < N; ++i) {
            d_[i] /= v;
        }

        return *this;
    }

    WPoint<T, N>& wmul(const_reference v) noexcept
    {
        for (size_type i = 0; i < N; ++i) {
            d_[i] *= v;
        }

        return *this;
    }

    Proj p() const __GM_NOEXCEPT_RELEASE__
    {
        auto& ww = w();
        check_ifd(!gm::cmp::zero(ww), "Zero weight");

        Proj result;
        for (size_type i = 0; i < N; ++i) {
            result[i] = d_[i] / ww;
        }
        return result;
    }

    const_reference w() const noexcept
    {
        return (*this)[N];
    }

    WPoint& operator+=(const WPoint& rhs) noexcept
    {
        for (size_type i = 0; i < N + 1; ++i) {
            d_[i] += rhs[i];
        }
        return *this;
    }

    WPoint& operator-=(const WPoint& rhs) noexcept
    {
        for (size_type i = 0; i < N + 1; ++i) {
            d_[i] -= rhs[i];
        }
        return *this;
    }

    WPoint& operator*=(const_reference rhs) noexcept
    {
        for (size_type i = 0; i < N + 1; ++i) {
            d_[i] *= rhs;
        }
        return *this;
    }

    WPoint& operator/=(const_reference rhs) __GM_NOEXCEPT_RELEASE__
    {
        check_ifd(!gm::cmp::zero(rhs), "Division by zero");

        for (size_type i = 0; i < N + 1; ++i) {
            d_[i] /= rhs;
        }
        return *this;
    }
};

template <class T, size_t N>
WPoint<T, 1> wdot(const WPoint<T, N>& lhs, const WPoint<T, N>& rhs) noexcept
{
    T wp {};
    for (size_t i = 0; i < N; ++i) {
        wp += lhs[i] * rhs[i];
    }
    return WPoint<T, 1>({wp, lhs[N] * rhs[N]});
}

template <class T, size_t N>
std::ostream& operator<<(std::ostream& os, const WPoint<T, N>& obj)
{
    fmt::print(os, "{{ \"wp\": {}, \"w\": {} }}",
               RangePrint(obj.begin(), std::prev(obj.end())), obj[N]);
    return os;
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
} // namespace gm

#endif // GEOM_MODEL_SRC_BSPLINE_WPOINT_HPP_