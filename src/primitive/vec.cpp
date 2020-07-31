#include <gm/point.hpp>
#include <gm/vec.hpp>

#include <gm/compare.hpp>
#include <util/math.hpp>

#include <fmt/ostream.h>

namespace gm {

Vec::Vec() noexcept
    : data_ {0, 0, 0}
{
}

Vec::Vec(value_type x, value_type y, value_type z) noexcept
    : data_ {x, y, z}
{
}

Vec::Vec(value_type magnitude, const std::array<value_type, N>& dir) noexcept
    : data_ {magnitude * dir[0], magnitude * dir[1], magnitude * dir[2]}
{
}

Vec::Vec(const std::array<value_type, N>& coord) noexcept
    : data_ {coord[0], coord[1], coord[2]}
{
}

Vec::Vec(const std::initializer_list<value_type>& list) noexcept
    : Vec()
{
    size_type i = 0;
    for (auto j = std::begin(list); i < N && j != std::end(list); ++j) {
        data_[i++] = *j;
    }
}

Vec::Vec(const Point& p) noexcept
    : data_ {p[0], p[1], p[2]}
{
}

Vec::Vec(const Point& lhs, const Point& rhs) noexcept
    : Vec(rhs - lhs)
{
}

Vec::pointer Vec::data() noexcept
{
    return data_;
}

Vec::const_pointer Vec::data() const noexcept
{
    return data_;
}

Vec::size_type Vec::size() const noexcept
{
    return N;
}

Vec::reference Vec::operator[](size_type i) noexcept
{
    return data_[i];
}

Vec::const_reference Vec::operator[](size_type i) const noexcept
{
    return data_[i];
}

Vec::iterator Vec::begin() noexcept
{
    return std::begin(data_);
}

Vec::iterator Vec::end() noexcept
{
    return std::end(data_);
}

Vec::const_iterator Vec::begin() const noexcept
{
    return std::begin(data_);
}

Vec::const_iterator Vec::end() const noexcept
{
    return std::end(data_);
}

std::array<Vec::value_type, Vec::N> Vec::raw() const noexcept
{
    return {data_[0], data_[1], data_[2]};
}

Vec& Vec::operator+=(const Vec& rhs) noexcept
{
    data_[0] += rhs[0];
    data_[1] += rhs[1];
    data_[2] += rhs[2];
    return *this;
}

Vec& Vec::operator-=(const Vec& rhs) noexcept
{
    data_[0] -= rhs[0];
    data_[1] -= rhs[1];
    data_[2] -= rhs[2];
    return *this;
}

Vec& Vec::operator*=(const_reference rhs) noexcept
{
    data_[0] *= rhs;
    data_[1] *= rhs;
    data_[2] *= rhs;
    return *this;
}

Vec& Vec::operator/=(const_reference rhs) __GM_NOEXCEPT_RELEASE__
{
    check_ifd(!cmp::zero(rhs), "Division by zero");

    data_[0] /= rhs;
    data_[1] /= rhs;
    data_[2] /= rhs;
    return *this;
}

Vec operator-(const Vec& obj) noexcept
{
    return {-obj[0], -obj[1], -obj[2]};
}

Vec operator+(const Vec& lhs, const Vec& rhs) noexcept
{
    auto result = lhs;
    return (result += rhs);
}

Vec operator-(const Vec& lhs, const Vec& rhs) noexcept
{
    auto result = lhs;
    return (result -= rhs);
}

Vec operator*(const Vec& lhs, Vec::const_reference rhs) noexcept
{
    auto result = lhs;
    return (result *= rhs);
}

Vec operator*(Vec::const_reference lhs, const Vec& rhs) noexcept
{
    auto result = rhs;
    return (result *= lhs);
}

Vec operator/(const Vec& lhs, Vec::const_reference rhs) __GM_NOEXCEPT_RELEASE__
{
    auto result = lhs;
    return (result /= rhs);
}

bool operator==(const Vec& lhs, const Vec& rhs) noexcept
{
    return cmp::near(lhs, rhs);
}

bool operator!=(const Vec& lhs, const Vec& rhs) noexcept
{
    return !(lhs == rhs);
}

Vec cross(const Vec& a, const Vec& b) noexcept
{
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]};
}

double sqr(const Vec& obj) noexcept
{
    return dot(obj, obj);
}

double norm(const Vec& obj) noexcept
{
    return sqrt(sqr(obj));
}

double dist(const Vec& lhs, const Vec& rhs) noexcept
{
    return norm(rhs - lhs);
}

bool isnan(const Vec& obj) noexcept
{
    return std::isnan(obj[0]) || std::isnan(obj[1]) || std::isnan(obj[2]);
}

bool isinf(const Vec& obj) noexcept
{
    return std::isinf(obj[0]) || std::isinf(obj[1]) || std::isinf(obj[2]);
}

double cos(const Vec& a, const Vec& b) noexcept
{
    return dot(a, b) / (norm(a) * norm(b));
}

double sin(const Vec& a, const Vec& b) noexcept
{
    return norm(cross(a, b)) / (norm(a) * norm(b));
}

double angle(const Vec& a, const Vec& b) noexcept
{
    auto x = dot(a, b);
    auto y = norm(cross(a, b));
    return atan2(y, x);
}

double mixed(const Vec& a, const Vec& b, const Vec& c) noexcept
{
    return dot(a, cross(b, c));
}

Vec unit(const Vec& obj) __GM_NOEXCEPT_RELEASE__
{
    auto n = norm(obj);
    check_ifd(!cmp::zero(n), "Unit vector of zero");

    return obj / n;
}

Vec bisect(const Vec& lhs, const Vec& rhs) __GM_NOEXCEPT_RELEASE__
{
    return unit(unit(lhs) + unit(rhs));
}

Vec bisect(std::initializer_list<Vec> list) __GM_NOEXCEPT_RELEASE__
{
    auto i = std::begin(list);
    auto result = *i++;
    for (; i != std::end(list); ++i) {
        result = unit(result) + unit(*i);
    }
    return unit(result);
}

std::ostream& operator<<(std::ostream& os, const Vec& obj)
{
    fmt::print(os, "[{:.5g}, {:.5g}, {:.5g}]", obj[0], obj[1], obj[2]);
    return os;
}

namespace cmp {
    bool isnear(const Vec& lhs, const Vec& rhs, Tolerance tol) noexcept
    {
        return zero(dist(lhs, rhs), tol);
    }

    bool near(const Vec& lhs, const Vec& rhs, Tolerance tol) noexcept
    {
        return isnear(lhs, rhs, tol);
    }

    bool zero(const Vec& lhs, Tolerance tol) noexcept
    {
        return zero(norm(lhs), tol);
    }
} // namespace cmp

} // namespace gm