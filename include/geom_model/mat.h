#ifndef GEOM_MODEL_INCLUDE_MAT_H_
#define GEOM_MODEL_INCLUDE_MAT_H_

#include "point.h"
#include "vec.h"

#include <memory>
#include <ostream>

class Mat;

Vec dot(const Mat& a, const Vec& v);
Point dot(const Mat& a, const Point& p);
Mat operator*(double x, const Mat& a);

std::ostream& operator<<(std::ostream& os, const Mat& mat);

class Mat {
    friend class Axis;

public:
    enum class RotateAx { X, Y, Z };
    enum class Direction : bool {
        CLOCKWISE = false,
        COUNTER_CLOCKWISE = true
    };

    static constexpr auto row_size = 3;
    static constexpr auto mat_size = 9;
    using data_type = std::array<double, mat_size>;

    ~Mat();
    Mat(Mat&& other) noexcept;
    Mat& operator=(Mat&& other) noexcept;
    Mat(const Mat& other);
    Mat& operator=(const Mat& other);

    Mat();
    explicit Mat(const data_type& coord);

    static Mat eye();
    static Mat rotate_x(double angle);
    static Mat rotate_y(double angle);
    static Mat rotate_z(double angle);
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
    static Mat rotate_x(double c, double s);
    static Mat rotate_y(double c, double s);
    static Mat rotate_z(double c, double s);

    struct Impl;
    std::unique_ptr<Impl> pimpl_;
};

namespace std {
template <>
struct hash<Mat> {
    size_t operator()(const Mat& key) const
    {
        return hasher_(key.raw());
    }

private:
    hash<Mat::data_type> hasher_;
};
}

#endif // GEOM_MODEL_INCLUDE_MAT_H_
