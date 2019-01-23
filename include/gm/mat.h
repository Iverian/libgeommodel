#ifndef GEOM_MODEL_INCLUDE_MAT_H_
#define GEOM_MODEL_INCLUDE_MAT_H_

#include "point.h"
#include "vec.h"

#include <memory>
#include <ostream>

namespace gm {

class Mat {
    friend class Axis;

public:
    static constexpr auto row_size = Vec::N;
    static constexpr auto mat_size = row_size * row_size;
    using data_type = std::array<double, mat_size>;

    ~Mat();
    Mat(Mat&& other) noexcept;
    Mat& operator=(Mat&& other) noexcept;
    Mat(const Mat& other);
    Mat& operator=(const Mat& other);

    Mat();
    explicit Mat(const data_type& coord);

    static Mat eye();
    static Mat rotate(double angle, const Vec& ax);

    const double& operator()(size_t i, size_t j) const noexcept;
    const double& operator[](size_t i) const noexcept;
    size_t size() const noexcept;
    const data_type& raw() const noexcept;

    Mat& operator+=(const Mat& other);
    Mat& operator-=(const Mat& other);
    Mat& operator*=(const Mat& other);
    Mat& operator*=(double x);
    Mat& operator/=(double x);
    Mat operator-() const;
    Mat operator+(const Mat& other) const;
    Mat operator-(const Mat& other) const;
    Mat operator*(double x) const;

    Vec operator*(const Vec& v) const;
    Point operator*(const Point& p) const;
    Mat operator*(const Mat& other) const;

    friend bool operator==(const Mat& lhs, const Mat& rhs);
    friend bool operator!=(const Mat& lhs, const Mat& rhs);

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl_;
};

Vec dot(const Mat& a, const Vec& v);
Point dot(const Mat& a, const Point& p);
Mat operator*(double x, const Mat& a);
std::ostream& operator<<(std::ostream& os, const Mat& mat);

} // namespace gm

namespace std {
template <>
struct hash<gm::Mat> {
    size_t operator()(const gm::Mat& key) const
    {
        size_t result = 0;
        for (size_t i = 0; i < gm::Mat::mat_size; ++i) {
            result = (result << 1) ^ hasher_(key[i]);
        }
        return result;
    }

private:
    hash<double> hasher_;
};
}

#endif // GEOM_MODEL_INCLUDE_MAT_H_
