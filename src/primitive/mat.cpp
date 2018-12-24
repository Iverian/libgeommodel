#include <geom_model/geom_util.h>
#include <geom_model/mat.h>
#include <util/math.h>

using namespace std;

struct Mat::Impl {
    Impl();
    explicit Impl(const Mat::data_type& coord);

    const double& call(size_t i, size_t j) const noexcept;
    const double& index(size_t i) const noexcept;
    size_t size() const noexcept;
    const data_type& raw() const noexcept;

    Impl& operator+=(const Impl& other);
    Impl& operator-=(const Impl& other);
    Impl& operator*=(const Impl& other);
    Impl& operator*=(double x);
    Impl& operator/=(double x);
    Impl operator-() const;
    Impl operator+(const Impl& other) const;
    Impl operator-(const Impl& other) const;
    Impl operator*(double x) const;

    array<double, 3> dot(const array<double, 3>& x);

    bool operator==(const Impl& other) const;

private:
    Mat::data_type coord_;
};

Mat::~Mat() = default;
Mat::Mat(Mat&& other) noexcept = default;
Mat& Mat::operator=(Mat&& other) noexcept = default;
Mat::Mat(const Mat& other)
    : pimpl_(make_unique<Impl>(*other.pimpl_))
{
}
Mat& Mat::operator=(const Mat& other)
{
    *pimpl_ = *other.pimpl_;
    return *this;
}

Mat::Mat()
    : pimpl_(make_unique<Mat::Impl>())
{
}

Mat::Mat(const Mat::data_type& coord)
    : pimpl_(make_unique<Mat::Impl>(coord))
{
}

Mat Mat::eye()
{
    return Mat({1, 0, 0, 0, 1, 0, 0, 0, 1});
}

Mat Mat::rotate_x(double angle)
{
    return rotate_x(cos(angle), sin(angle));
}

Mat Mat::rotate_y(double angle)
{
    return rotate_y(cos(angle), sin(angle));
}

Mat Mat::rotate_z(double angle)
{
    return rotate_z(cos(angle), sin(angle));
}

Mat Mat::rotate_x(double c, double s)
{
    return Mat({1, 0, 0, 0, c, s, 0, -s, c});
}

Mat Mat::rotate_y(double c, double s)
{
    return Mat({c, 0, s, 0, 1, 0, -s, 0, c});
}

Mat Mat::rotate_z(double c, double s)
{
    return Mat({c, s, 0, -s, c, 0, 0, 0, 1});
}

#define sqr(x) ((x) * (x))

Mat Mat::rotate(double angle, const Vec& ax)
{
    auto u = unit(ax);
    auto x = u[0], y = u[1], z = u[2];
    auto c = cos(angle), s = sin(angle);
    return Mat({c + sqr(x) * (1 - c), x * y * (1 - c) - z * s,
                x * z * (1 - c) + y * s, y * x * (1 - c) + z * s,
                c + sqr(y) * (1 - c), y * z * (1 - c) - x * s,
                z * x * (1 - c) - y * s, z * y * (1 - c) + x * s,
                c + sqr(z) * (1 - c)});
}

#undef sqr

const double& Mat::operator()(size_t i, size_t j) const noexcept
{
    return pimpl_->call(i, j);
}

const double& Mat::operator[](size_t i) const noexcept
{
    return pimpl_->index(i);
}

size_t Mat::size() const noexcept
{
    return pimpl_->size();
}

const Mat::data_type& Mat::raw() const noexcept
{
    return pimpl_->raw();
}

Mat& Mat::operator+=(const Mat& other)
{
    *pimpl_ += *other.pimpl_;
    return *this;
}

Mat& Mat::operator-=(const Mat& other)
{
    *pimpl_ -= *other.pimpl_;
    return *this;
}

Mat& Mat::operator*=(const Mat& other)
{
    *pimpl_ *= *other.pimpl_;
    return *this;
}

Mat& Mat::operator*=(double x)
{
    *pimpl_ *= x;
    return *this;
}

Mat& Mat::operator/=(double x)
{
    *pimpl_ /= x;
    return *this;
}

Mat Mat::operator-() const
{
    Mat result = *this;
    return (result *= -1);
}

Mat Mat::operator+(const Mat& other) const
{
    Mat result = *this;
    return (result += other);
}

Mat Mat::operator-(const Mat& other) const
{
    Mat result = *this;
    return (result -= other);
}

Mat Mat::operator*(double x) const
{
    Mat result = *this;
    return (result *= x);
}

Vec Mat::operator*(const Vec& v) const
{
    return Vec(pimpl_->dot(v.raw()));
}

Point Mat::operator*(const Point& p) const
{
    return Point(pimpl_->dot(p.raw()));
}

Mat Mat::operator*(const Mat& other) const
{
    auto result = *this;
    return (result *= other);
}

bool operator==(const Mat& lhs, const Mat& rhs)
{
    return *lhs.pimpl_ == *rhs.pimpl_;
}

bool operator!=(const Mat& lhs, const Mat& rhs)
{
    return !(lhs == rhs);
}

Vec dot(const Mat& a, const Vec& v)
{
    return a * v;
}

Point dot(const Mat& a, const Point& p)
{
    return a * p;
}

Mat operator*(double x, const Mat& a)
{
    return a * x;
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

Mat::Impl::Impl()
    : coord_{}
{
}

Mat::Impl::Impl(const Mat::data_type& coord)
    : coord_{coord}
{
}

const double& Mat::Impl::call(size_t i, size_t j) const noexcept
{
    return coord_[row_size * i + j];
}

const double& Mat::Impl::index(size_t i) const noexcept
{
    return coord_[i];
}

size_t Mat::Impl::size() const noexcept
{
    return mat_size;
}

Mat::Impl& Mat::Impl::operator+=(const Mat::Impl& other)
{
    for (size_t i = 0; i < mat_size; ++i)
        coord_[i] += other.coord_[i];
    return *this;
}

Mat::Impl& Mat::Impl::operator-=(const Mat::Impl& other)
{
    for (size_t i = 0; i < mat_size; ++i)
        coord_[i] -= other.coord_[i];
    return *this;
}

Mat::Impl& Mat::Impl::operator*=(double x)
{
    for (size_t i = 0; i < mat_size; ++i)
        coord_[i] *= x;
    return *this;
}

#define _(i, j) ((i) * (row_size) + (j))

Mat::Impl& Mat::Impl::operator*=(const Impl& other)
{
    Mat::data_type result;
    for (size_t i = 0; i < row_size; ++i) {
        for (size_t j = 0; j < row_size; ++j) {
            result[_(i, j)] = 0;
            for (size_t k = 0; k < row_size; ++k) {
                result[_(i, j)] += coord_[_(i, k)] * other.coord_[_(k, j)];
            }
        }
    }
    coord_ = result;
    return *this;
}

#undef _

Mat::Impl& Mat::Impl::operator/=(double x)
{
    for (size_t i = 0; i < mat_size; ++i)
        coord_[i] /= x;
    return *this;
}

Mat::Impl Mat::Impl::operator-() const
{
    Mat::Impl result = *this;
    return (result *= -1);
}

Mat::Impl Mat::Impl::operator+(const Mat::Impl& other) const
{
    Mat::Impl result = *this;
    return (result += other);
}

Mat::Impl Mat::Impl::operator-(const Mat::Impl& other) const
{
    Mat::Impl result = *this;
    return (result -= other);
}

Mat::Impl Mat::Impl::operator*(double x) const
{
    Mat::Impl result = *this;
    return result;
}

array<double, 3> Mat::Impl::dot(const array<double, 3>& x)
{
    array<double, 3> result{};
    for (size_t i = 0; i < row_size; ++i) {
        for (size_t j = 0; j < row_size; ++j) {
            result[i] += call(i, j) * x[j];
        }
    }
    return result;
}

const Mat::data_type& Mat::Impl::raw() const noexcept
{
    return coord_;
}

bool Mat::Impl::operator==(const Mat::Impl& other) const
{
    return coord_ == other.coord_;
}
