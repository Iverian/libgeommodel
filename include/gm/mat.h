#ifndef GEOM_MODEL_INCLUDE_GM_MAT_H_
#define GEOM_MODEL_INCLUDE_GM_MAT_H_

#include "point.h"
#include "vec.h"

#include <memory>
#include <ostream>
#include <optional>

namespace gm {

class Mat {
    friend class Axis;

public:
    static constexpr auto row_size = Vec::N;
    static constexpr auto mat_size = row_size * row_size;
    using data_type = std::array<double, mat_size>;

    Mat();
    explicit Mat(const data_type& coord);

    static Mat eye();
    static Mat rotate(double angle, const Vec& ax);

    double& operator()(size_t i, size_t j) noexcept;
    double& operator[](size_t i) noexcept;
    const double& operator()(size_t i, size_t j) const noexcept;
    const double& operator[](size_t i) const noexcept;
    size_t size() const noexcept;
    const data_type& raw() const noexcept;

    Mat& operator+=(const Mat& rhs);
    Mat& operator-=(const Mat& rhs);
    Mat& operator*=(const Mat& rhs);
    Mat& operator*=(double rhs);
    Mat& operator/=(double rhs);

    friend Mat operator+(const Mat& lhs, const Mat& rhs);
    friend Mat operator-(const Mat& lhs, const Mat& rhs);
    friend Mat operator*(const Mat& lhs, double rhs);
    friend Mat operator*(double lhs, const Mat& rhs);
    friend Mat operator/(const Mat& lhs, double rhs);

    friend bool operator==(const Mat& lhs, const Mat& rhs);
    friend bool operator!=(const Mat& lhs, const Mat& rhs);

private:
    data_type data_;
};

double det(const Mat& x) noexcept;
std::optional<Mat> inverse(const Mat& x) noexcept;
Vec dot(const Mat& a, const Vec& v);
Point dot(const Mat& a, const Point& p);
std::ostream& operator<<(std::ostream& os, const Mat& mat);

} // namespace gm

namespace std {
template <>
struct hash<gm::Mat> {
    size_t operator()(const gm::Mat& key) const
    {
        size_t result = 16651;
        for (size_t i = 0; i < gm::Mat::mat_size; ++i) {
            result = (result << 1) ^ hasher_(key[i]);
        }
        return result;
    }

private:
    hash<double> hasher_;
};
}

#endif // GEOM_MODEL_INCLUDE_GM_MAT_H_
