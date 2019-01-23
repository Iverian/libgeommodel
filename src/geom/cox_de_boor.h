#ifndef GEOM_MODEL_SRC_GEOM_COX_DE_BOOR_H_
#define GEOM_MODEL_SRC_GEOM_COX_DE_BOOR_H_

#include <gm/point.h>
#include <util/debug.h>
#include <util/math.h>

#include <array>
#include <type_traits>
#include <vector>

template <class T>
struct VectorView {
    using value_type = typename std::vector<T>::value_type;
    using size_type = typename std::vector<T>::size_type;
    using const_reference = typename std::vector<T>::const_reference;
    using const_pointer = typename std::vector<T>::const_pointer;

    const_reference operator[](size_type i) const noexcept;
    const_reference at(size_type i) const;
    const_reference front() const;
    const_reference back() const;
    size_type size() const;

    VectorView();
    VectorView(const std::vector<T>& source);
    VectorView(const_pointer data, size_type size);

    const_pointer data_;
    size_type size_;
};


template <class T>
struct CoxDeBoor {
    using value_type = typename std::decay<T>::type;
    using reference = value_type&;
    using const_reference = const value_type&;

    friend struct Proxy;

    struct Proxy {
        friend struct CoxDeBoor;

        value_type get(size_t k) const;
        std::vector<value_type> range(size_t k) const;
        Proxy(const CoxDeBoor& parent, double t, size_t p);

    private:
        const CoxDeBoor* parent_;
        double t_;
        size_t p_;
    };

    CoxDeBoor();
    CoxDeBoor(size_t order, const std::vector<double>& knots,
              const std::vector<value_type>& cpoints);
    CoxDeBoor(size_t order, VectorView<double> knots,
              VectorView<value_type> cpoints);

    Proxy proxy(double t) const;
    double pad(double t) const;
    size_t interval(double t) const;
    value_type pget(size_t i, size_t k) const;
    value_type eval(double t, size_t p, std::vector<value_type>& cp) const;

private:
    size_t order_;
    VectorView<double> knots_;
    VectorView<value_type> cpoints_;
};

template <class T>
CoxDeBoor<T>::CoxDeBoor()
    : order_(0)
    , knots_()
    , cpoints_()
{
}

template <class T>
CoxDeBoor<T>::CoxDeBoor(size_t order, const std::vector<double>& knots,
                        const std::vector<value_type>& cpoints)
    : order_(order)
    , knots_(knots)
    , cpoints_(cpoints)
{
}

template <class T>
CoxDeBoor<T>::CoxDeBoor(size_t order, VectorView<double> knots,
                        VectorView<value_type> cpoints)
    : order_(std::move(order))
    , knots_(std::move(knots))
    , cpoints_(std::move(cpoints))
{
}

template <class T>
typename CoxDeBoor<T>::Proxy CoxDeBoor<T>::proxy(double t) const
{
    t = pad(t);
    return Proxy(*this, t, interval(t));
}

template <class T>
double CoxDeBoor<T>::pad(double t) const
{
    return ::pad(t, knots_.front(), knots_.back());
}

template <class T>
size_t CoxDeBoor<T>::interval(double t) const
{
    size_t result = 0;
    size_t end = knots_.size() - 1;
    for (size_t i = 0; i < end; ++i) {
        if (knots_[i] <= t && t < knots_[i + 1]) {
            result = i;
            break;
        }
    }
    return result;
}

template <class T>
typename CoxDeBoor<T>::value_type CoxDeBoor<T>::pget(size_t i, size_t k) const
{
    value_type result;
    if (k == 0) {
        result = cpoints_.at(i);
    } else {
        result = double(order_ - k) * (pget(i, k - 1) - pget(i - 1, k - 1))
            / (knots_.at(i + order_ - 1) - knots_.at(i));
    }
    return result;
}

template <class T>
typename CoxDeBoor<T>::value_type
CoxDeBoor<T>::eval(double t, size_t p, std::vector<value_type>& cp) const
{
    auto n = cp.size();
    for (size_t i = 1; i < n; ++i) {
        for (size_t j = n - 1; j >= i; --j) {
            auto k = p - n + 1 + j;
            auto a = t - knots_[k];
            auto b = knots_[k + n - i] - t;
            cp[j] = (cp[j] * a + cp[j - 1] * b) / (a + b);
        }
    }
    return cp.back();
}

template <class T>
typename CoxDeBoor<T>::value_type CoxDeBoor<T>::Proxy::get(size_t k) const
{
    value_type result;
    if (k >= parent_->order_) {
        result = value_type();
    } else {
        auto n = parent_->order_ - k;
        std::vector<value_type> cp(n);

        for (size_t i = 0; i < n; ++i) {
            cp[i] = parent_->pget(p_ - n + 1 + i, k);
        }
        result = parent_->eval(t_, p_, cp);
    }

    return result;
}

template <class T>
std::vector<typename CoxDeBoor<T>::value_type>
CoxDeBoor<T>::Proxy::range(size_t k) const
{
    std::vector<value_type> result(k);
    for (size_t i = 0; i < k; ++i) {
        result[i] = get(i);
    }
    return result;
}

template <class T>
CoxDeBoor<T>::Proxy::Proxy(const CoxDeBoor<T>& parent, double t, size_t p)
    : parent_(&parent)
    , t_(t)
    , p_(p)
{
}

template <class T>
VectorView<T>::VectorView()
    : data_(nullptr)
    , size_(0)
{
}

template <class T>
VectorView<T>::VectorView(const std::vector<T>& source)
    : data_(source.data())
    , size_(source.size())
{
}

template <class T>
VectorView<T>::VectorView(const_pointer data, size_type size)
    : data_(std::move(data))
    , size_(std::move(size))
{
}

template <class T>
typename VectorView<T>::const_reference VectorView<T>::
operator[](size_type i) const noexcept
{
    return data_[i];
}

template <class T>
typename VectorView<T>::const_reference VectorView<T>::at(size_type i) const
{
    CHECK_(i < size_, std::out_of_range, "i: {} is out of range [0;{})", i,
           size_);
    return data_[i];
}

template <class T>
typename VectorView<T>::const_reference VectorView<T>::front() const
{
    return data_[0];
}

template <class T>
typename VectorView<T>::const_reference VectorView<T>::back() const
{
    return data_[size_ - 1];
}

template <class T>
typename VectorView<T>::size_type VectorView<T>::size() const
{
    return size_;
}

#endif // GEOM_MODEL_SRC_GEOM_COX_DE_BOOR_H_