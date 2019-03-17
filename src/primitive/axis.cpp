#include <gm/axis.hpp>
#include <gm/mat.hpp>

#include <fmt/ostream.h>

namespace gm {

Axis::Axis() noexcept
    : c_(0, 0, 0)
    , basis_ {Vec(1, 0, 0), Vec(0, 1, 0), Vec(0, 0, 1)}
{
}

Axis::Axis(const Point& c, const std::array<Vec, N>& basis) noexcept
    : c_(c)
    , basis_ {basis[0], basis[1], basis[2]}
{
    std::array<double, 2> coeff;
    for (size_t i = 1; i < 3; ++i) {
        coeff[i - 1] = -dot(basis_[i - 1], basis_[i]) / sqr(basis_[i - 1]);
        for (size_t j = 0; j < i; ++j) {
            basis_[i] += coeff[j] * basis_[j];
        }
    }
    for (auto& i : basis_)
        i = unit(i);
}

const Vec& Axis::operator[](size_t i) const noexcept
{
    return basis_[i];
}

const Vec& Axis::x() const noexcept
{
    return basis_[0];
}

const Vec& Axis::y() const noexcept
{
    return basis_[1];
}

const Vec& Axis::z() const noexcept
{
    return basis_[2];
}

const Point& Axis::c() const noexcept
{
    return c_;
}

Axis::View Axis::get_view() const noexcept
{
    return std::make_tuple(c_, basis_[0], basis_[1], basis_[2]);
}

Vec Axis::rotate_x(double angle, const Vec& v) const noexcept
{
    return dot(Mat::rotate(angle, basis_[0]), v);
}

Vec Axis::rotate_y(double angle, const Vec& v) const noexcept
{
    return dot(Mat::rotate(angle, basis_[1]), v);
}

Vec Axis::rotate_z(double angle, const Vec& v) const noexcept
{
    return dot(Mat::rotate(angle, basis_[2]), v);
}

Point Axis::rotate_x(double angle, const Point& p) const noexcept
{
    return dot(Mat::rotate(angle, basis_[0]), p - c_) + c_;
}

Point Axis::rotate_y(double angle, const Point& p) const noexcept
{
    return dot(Mat::rotate(angle, basis_[1]), p - c_) + c_;
}

Point Axis::rotate_z(double angle, const Point& p) const noexcept
{
    return dot(Mat::rotate(angle, basis_[2]), p - c_) + c_;
}

Point Axis::global(const Point& p) const noexcept
{
    return p[0] * basis_[0] + p[1] * basis_[1] + p[2] * basis_[2] + c_;
}

Vec Axis::global(const Vec& v) const noexcept
{
    return v[0] * basis_[0] + v[1] * basis_[1] + v[2] * basis_[2];
}

Point Axis::pglobal(double x, double y, double z) const noexcept
{
    return global(Point(x, y, z));
}

Vec Axis::vglobal(double x, double y, double z) const noexcept
{
    return global(Vec(x, y, z));
}

Axis Axis::from_xy(const Vec& x, const Vec& ref_y,
                   const Point& center) noexcept
{
    auto y = ref_y - dot(ref_y, x) * x;
    auto z = cross(x, y);
    return Axis(center, {x, y, z});
}

Axis Axis::from_zx(const Vec& z, const Vec& ref_x,
                   const Point& center) noexcept
{
    auto x = ref_x - dot(ref_x, z) * z;
    auto y = cross(z, x);
    return Axis(center, {x, y, z});
}

Axis Axis::from_abc(const Point& a, const Point& b, const Point& c) noexcept
{
    return Axis::from_xy(Vec(a - b), Vec(c - b), b);
}

bool operator==(const Axis& lhs, const Axis& rhs) noexcept
{
    return (lhs.c_ == rhs.c_) && (lhs[0] == rhs[0]) && (lhs[1] == rhs[1])
        && (lhs[2] == rhs[2]);
}

bool operator!=(const Axis& lhs, const Axis& rhs) noexcept
{
    return !(lhs == rhs);
}

std::ostream& operator<<(std::ostream& os, const Axis& ax)
{
    fmt::print(os, "{{ \"x\": {0}, \"y\": {1}, \"z\": {2}, \"c\": {3} }}",
               ax[0], ax[1], ax[2], ax.c());
    return os;
}

} // namespace gm