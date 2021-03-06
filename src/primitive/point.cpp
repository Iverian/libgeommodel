#include <gm/point.hpp>
#include <gm/vec.hpp>

#include <fmt/ostream.h>
#include <util/debug.hpp>
#include <util/math.hpp>

namespace gm {

Point::Point() noexcept
    : data_ {0, 0, 0}
{
}

Point::Point(value_type x, value_type y, value_type z) noexcept
    : data_ {x, y, z}
{
}

Point::Point(const std::array<value_type, N>& coord) noexcept
    : data_ {coord[0], coord[1], coord[2]}
{
}

Point::Point(const std::initializer_list<value_type>& list) noexcept
    : Point()
{
    size_type i = 0;
    for (auto j = std::begin(list); i < N && j != std::end(list); ++i, ++j) {
        data_[i] = *j;
    }
}

Point::Point(const Vec& lhs, const Vec& rhs) noexcept
    : Point(rhs - lhs)
{
}

Point::Point(const Vec& v) noexcept
    : data_ {v[0], v[1], v[2]}
{
}

Point::pointer Point::data() noexcept
{
    return data_;
}

Point::const_pointer Point::data() const noexcept
{
    return data_;
}

Point::size_type Point::size() const noexcept
{
    return N;
}

Point::reference Point::operator[](size_type i) noexcept
{
    return data_[i];
}

Point::const_reference Point::operator[](size_type i) const noexcept
{
    return data_[i];
}

Point::iterator Point::begin() noexcept
{
    return std::begin(data_);
}

Point::iterator Point::end() noexcept
{
    return std::end(data_);
}

Point::const_iterator Point::begin() const noexcept
{
    return std::begin(data_);
}

Point::const_iterator Point::end() const noexcept
{
    return std::end(data_);
}

std::array<Point::value_type, Point::N> Point::raw() const noexcept
{
    return {data_[0], data_[1], data_[2]};
}

Point& Point::operator+=(const Point& rhs) noexcept
{
    data_[0] += rhs[0];
    data_[1] += rhs[1];
    data_[2] += rhs[2];
    return *this;
}

Point& Point::operator-=(const Point& rhs) noexcept
{
    data_[0] -= rhs[0];
    data_[1] -= rhs[1];
    data_[2] -= rhs[2];
    return *this;
}

Point& Point::operator*=(const_reference rhs) noexcept
{
    data_[0] *= rhs;
    data_[1] *= rhs;
    data_[2] *= rhs;
    return *this;
}

Point& Point::operator/=(const_reference rhs) __GM_NOEXCEPT_RELEASE__
{
    check_ifd(!cmp::zero(rhs), "Division by zero");

    data_[0] /= rhs;
    data_[1] /= rhs;
    data_[2] /= rhs;
    return *this;
}

Point operator+(const Point& lhs, const Point& rhs) noexcept
{
    auto result = lhs;
    return (result += rhs);
}

Point operator-(const Point& lhs, const Point& rhs) noexcept
{
    auto result = lhs;
    return (result -= rhs);
}

Point operator*(const Point& lhs, Point::const_reference rhs) noexcept
{
    auto result = lhs;
    return (result *= rhs);
}

Point operator*(Point::const_reference lhs, const Point& rhs) noexcept
{
    auto result = rhs;
    return (result *= lhs);
}

Point operator/(const Point& lhs,
                Point::const_reference rhs) __GM_NOEXCEPT_RELEASE__
{
    auto result = lhs;
    return (result /= rhs);
}

Point operator-(const Point& obj) noexcept
{
    return {-obj[0], -obj[1], -obj[2]};
}

Point operator+(const Point& lhs, const Vec& rhs) noexcept
{
    auto result = lhs;
    result[0] += rhs[0];
    result[1] += rhs[1];
    result[2] += rhs[2];
    return result;
}

Point operator+(const Vec& lhs, const Point& rhs) noexcept
{
    auto result = rhs;
    result[0] += lhs[0];
    result[1] += lhs[1];
    result[2] += lhs[2];
    return result;
}

Point operator-(const Point& lhs, const Vec& rhs) noexcept
{
    auto result = lhs;
    result[0] -= rhs[0];
    result[1] -= rhs[1];
    result[2] -= rhs[2];
    return result;
}

Point operator-(const Vec& lhs, const Point& rhs) noexcept
{
    auto result = lhs;
    result[0] -= rhs[0];
    result[1] -= rhs[1];
    result[2] -= rhs[2];
    return Point(result);
}

bool operator==(const Point& lhs, const Point& rhs) noexcept
{
    return cmp::near(lhs, rhs);
}

bool operator!=(const Point& lhs, const Point& rhs) noexcept
{
    return !(lhs == rhs);
}

double sqr(const Point& obj) noexcept
{
    return dot(obj, obj);
}

double dist(const Point& lhs, const Point& rhs) noexcept
{
    return sqrt(sqr(rhs - lhs));
}

bool isnan(const Point& obj) noexcept
{
    return std::isnan(obj[0]) || std::isnan(obj[1]) || std::isnan(obj[2]);
}

bool isinf(const Point& obj) noexcept
{
    return std::isinf(obj[0]) || std::isinf(obj[1]) || std::isinf(obj[2]);
}

std::ostream& operator<<(std::ostream& os, const Point& obj)
{
    fmt::print(os, "[{:.5g}, {:.5g}, {:.5g}]", obj[0], obj[1], obj[2]);
    return os;
}

namespace cmp {
    bool isnear(const Point& lhs, const Point& rhs, Tolerance tol) noexcept
    {
        return zero(dist(lhs, rhs), tol);
    }

    bool near(const Point& lhs, const Point& rhs, Tolerance tol) noexcept
    {
        return isnear(lhs, rhs, tol);
    }

    bool zero(const Point& lhs, Tolerance tol) noexcept
    {
        return zero(::sqrt(sqr(lhs)), tol);
    }
} // namespace cmp

} // namespace gm