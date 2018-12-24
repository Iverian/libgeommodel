#ifndef GEOM_MODEL_SRC_PRIMITIVE_VECIMPL_H_
#define GEOM_MODEL_SRC_PRIMITIVE_VECIMPL_H_

#include <array>
#include <ctime>
#include <functional>
#include <ostream>

class VecImpl;

double dot(const VecImpl& a, const VecImpl& b);
double sqr(const VecImpl& a);
std::ostream& operator<<(std::ostream& os, const VecImpl& v);

class VecImpl {
    using data_type = std::array<double, 3>;

public:
    VecImpl();
    VecImpl(double x, double y, double z);
    VecImpl(std::initializer_list<double> data);
    explicit VecImpl(const data_type& data);

    const double& operator[](size_t i) const noexcept;
    const double* data() const noexcept;
    size_t size() const noexcept;
    const data_type& raw() const noexcept;

    data_type::iterator begin();
    data_type::iterator end();
    data_type::const_iterator cbegin() const;
    data_type::const_iterator cend() const;

    VecImpl& operator+=(const VecImpl& other);
    VecImpl& operator-=(const VecImpl& other);
    VecImpl& operator*=(double x);
    VecImpl& operator/=(double x);

    double dot(const VecImpl& other) const;
    double sqr() const;
    bool isnan() const;

    friend bool operator==(const VecImpl& lhs, const VecImpl& rhs);
    friend bool operator!=(const VecImpl& lhs, const VecImpl& rhs);

private:
    data_type data_;
};

#endif // GEOM_MODEL_SRC_PRIMITIVE_VECIMPL_H_
