#ifndef GEOM_MODEL_INCLUDE_PARAMETRIC_POINT_H_
#define GEOM_MODEL_INCLUDE_PARAMETRIC_POINT_H_

#include <algorithm>
#include <iostream>

struct ParametricPoint;

std::ostream& operator<<(std::ostream& os, const ParametricPoint& p);

struct ParametricPoint {
    constexpr ParametricPoint(double pu, double pv)
        : u{pu}
        , v{pv}
    {
    }

    constexpr ParametricPoint(const ParametricPoint& other)
        : u{other.u}
        , v{other.v}
    {
    }

    constexpr ParametricPoint(ParametricPoint&& other) noexcept
        : u{other.u}
        , v{other.v}
    {
    }

    constexpr ParametricPoint& operator=(const ParametricPoint& other)
        = default;

    constexpr ParametricPoint& operator=(ParametricPoint&& other) noexcept
        = default;

    ~ParametricPoint() = default;

    constexpr ParametricPoint& operator+=(const ParametricPoint& other)
    {
        u += other.u;
        v += other.v;
        return *this;
    }

    constexpr ParametricPoint& operator-=(const ParametricPoint& other)
    {
        u -= other.u;
        v -= other.v;
        return *this;
    }

    constexpr ParametricPoint& operator*=(double x)
    {
        u *= x;
        v *= x;
        return *this;
    }

    constexpr ParametricPoint& operator/=(double x)
    {
        u /= x;
        v /= x;
        return *this;
    }

    constexpr ParametricPoint operator+(const ParametricPoint& other) const
    {
        auto result = *this;
        return (result += other);
    }
    constexpr ParametricPoint operator-(const ParametricPoint& other) const
    {
        auto result = *this;
        return (result -= other);
    }

    constexpr ParametricPoint operator*(double x) const
    {
        auto result = *this;
        return (result *= x);
    }

    constexpr ParametricPoint operator/(double x) const
    {
        auto result = *this;
        return (result /= x);
    }

    constexpr ParametricPoint()
        : u{0}
        , v{0}
    {
    }

    double u;
    double v;
};

constexpr ParametricPoint operator*(double x, const ParametricPoint& p)
{
    return (p * x);
}

constexpr bool operator==(const ParametricPoint& lhs,
                          const ParametricPoint& rhs)
{
    return (lhs.u == lhs.u) && (rhs.v == rhs.v);
}

constexpr bool operator!=(const ParametricPoint& lhs,
                          const ParametricPoint& rhs)
{
    return !(lhs == rhs);
}

constexpr ParametricPoint abs(const ParametricPoint& p)
{
    return {std::abs(p.u), std::abs(p.v)};
}

constexpr double max(const ParametricPoint& p)
{
    return std::max(p.u, p.v);
}

constexpr double min(const ParametricPoint& p)
{
    return std::min(p.u, p.v);
}

#endif // GEOM_MODEL_INCLUDE_PARAMETRIC_POINT_H_