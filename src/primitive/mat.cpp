#include <gm/compare.h>
#include <gm/mat.h>
#include <util/math.h>

#include <cmath>

using namespace std;

namespace gm {

#define _(i, j) ((i) * (row_size) + (j))

Mat Mat::eye()
{
    return Mat({1, 0, 0, 0, 1, 0, 0, 0, 1});
}

#define sqr(x) ((x) * (x))

Mat Mat::rotate(double angle, const Vec& ax)
{
    auto u = unit(ax);
    auto x = u[0], y = u[1], z = u[2];
    auto c = ::cos(angle), s = ::sin(angle);
    return Mat({c + sqr(x) * (1 - c), x * y * (1 - c) - z * s,
                x * z * (1 - c) + y * s, y * x * (1 - c) + z * s,
                c + sqr(y) * (1 - c), y * z * (1 - c) - x * s,
                z * x * (1 - c) - y * s, z * y * (1 - c) + x * s,
                c + sqr(z) * (1 - c)});
}

#undef sqr

double det(const Mat& x) noexcept
{
    return x(0, 0) * (x(1, 1) * x(2, 2) - x(1, 2) * x(2, 1))
        - x(0, 1) * (x(1, 0) * x(2, 2) - x(1, 2) * x(2, 0))
        + x(0, 2) * (x(1, 0) * x(2, 1) - x(1, 1) * x(2, 0));
}

::optional<Mat> inverse(const Mat& x) noexcept
{
    if (auto d = det(x); !iszero(d)) {
        Mat result;
        result(0, 0) = x(1, 1) * x(2, 2) - x(1, 2) * x(2, 1);
        result(0, 1) = x(0, 2) * x(2, 1) - x(0, 1) * x(2, 2);
        result(0, 2) = x(0, 1) * x(1, 2) - x(0, 2) * x(1, 1);
        result(1, 0) = x(1, 2) * x(2, 0) - x(1, 0) * x(2, 2);
        result(1, 1) = x(0, 0) * x(2, 2) - x(0, 2) * x(2, 0);
        result(1, 2) = x(0, 2) * x(1, 0) - x(0, 0) * x(1, 2);
        result(2, 0) = x(1, 0) * x(2, 1) - x(1, 1) * x(2, 0);
        result(2, 1) = x(0, 1) * x(2, 0) - x(0, 0) * x(2, 1);
        result(2, 2) = x(0, 0) * x(1, 1) - x(0, 1) * x(1, 0);

        return result / d;
    }
    return nullopt;
}

Vec dot(const Mat& a, const Vec& v)
{
    Vec result;
    result[0] = a(0, 0) * v[0] + a(0, 1) * v[1] + a(0, 2) * v[2];
    result[1] = a(1, 0) * v[0] + a(1, 1) * v[1] + a(1, 2) * v[2];
    result[2] = a(2, 0) * v[0] + a(2, 1) * v[1] + a(2, 2) * v[2];
    return result;
}

Point dot(const Mat& a, const Point& p)
{
    Point result;
    result[0] = a(0, 0) * p[0] + a(0, 1) * p[1] + a(0, 2) * p[2];
    result[1] = a(1, 0) * p[0] + a(1, 1) * p[1] + a(1, 2) * p[2];
    result[2] = a(2, 0) * p[0] + a(2, 1) * p[1] + a(2, 2) * p[2];
    return result;
}

ostream& operator<<(ostream& os, const Mat& mat)
{
    os << "[";
    for (size_t i = 0; i < mat.row_size; ++i) {
        os << "[";
        for (size_t j = 0; j < mat.row_size; ++j) {
            os << mat(i, j) << (j + 1 != mat.row_size ? ", " : "");
        }
        os << "]" << (i + 1 != mat.row_size ? ", " : "");
    }
    return os << "]";
}

Mat::Mat()
    : data_ {}
{
}

Mat::Mat(const Mat::data_type& data)
    : data_ {data}
{
}

double& Mat::operator()(size_t i, size_t j) noexcept
{
    return data_[_(i, j)];
}
double& Mat::operator[](size_t i) noexcept
{
    return data_[i];
}

const double& Mat::operator()(size_t i, size_t j) const noexcept
{
    return data_[_(i, j)];
}

const double& Mat::operator[](size_t i) const noexcept
{
    return data_[i];
}

size_t Mat::size() const noexcept
{
    return mat_size;
}

Mat& Mat::operator+=(const Mat& rhs)
{
    for (size_t i = 0; i < mat_size; ++i)
        data_[i] += rhs.data_[i];
    return *this;
}

Mat& Mat::operator-=(const Mat& other)
{
    for (size_t i = 0; i < mat_size; ++i)
        data_[i] -= other.data_[i];
    return *this;
}

Mat& Mat::operator*=(double x)
{
    for (size_t i = 0; i < mat_size; ++i)
        data_[i] *= x;
    return *this;
}

Mat& Mat::operator*=(const Mat& rhs)
{
    Mat::data_type result;
    for (size_t i = 0; i < row_size; ++i) {
        for (size_t j = 0; j < row_size; ++j) {
            result[_(i, j)] = 0;
            for (size_t k = 0; k < row_size; ++k) {
                result[_(i, j)] += data_[_(i, k)] * rhs.data_[_(k, j)];
            }
        }
    }
    data_ = result;
    return *this;
}

Mat& Mat::operator/=(double x)
{
    for (size_t i = 0; i < mat_size; ++i)
        data_[i] /= x;
    return *this;
}

const Mat::data_type& Mat::raw() const noexcept
{
    return data_;
}

Mat operator+(const Mat& lhs, const Mat& rhs)
{
    Mat result = lhs;
    return (result += rhs);
}

Mat operator-(const Mat& lhs, const Mat& rhs)
{
    Mat result = lhs;
    return (result -= rhs);
}

Mat operator*(const Mat& lhs, double rhs)
{
    Mat result = lhs;
    return (result *= rhs);
}

Mat operator*(double lhs, const Mat& rhs)
{
    Mat result = rhs;
    return (result *= lhs);
}

Mat operator/(const Mat& lhs, double rhs)
{
    Mat result = lhs;
    return (result /= rhs);
}

bool operator==(const Mat& lhs, const Mat& rhs)
{
    return lhs.data_ == rhs.data_;
}

bool operator!=(const Mat& lhs, const Mat& rhs)
{
    return !(lhs == rhs);
}

#undef _

} // namespace gm