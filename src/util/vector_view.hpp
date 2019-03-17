#ifndef GEOM_MODEL_SRC_UTIL_VECTOR_VIEW_HPP_
#define GEOM_MODEL_SRC_UTIL_VECTOR_VIEW_HPP_

#include "debug.hpp"

#include <vector>

template <class T, class Allocator = std::allocator<T>>
class VectorView {
    using owner_type = std::vector<T, Allocator>;

public:
    using value_type = typename owner_type::value_type;
    using size_type = typename owner_type::size_type;
    using const_reference = typename owner_type::const_reference;
    using const_pointer = typename owner_type::const_pointer;

    VectorView()
        : data_(nullptr)
        , size_(0)
    {
    }
    explicit VectorView(const owner_type& owner)
        : data_(owner.data())
        , size_(owner.size())
    {
    }
    VectorView(const_pointer data, size_type size)
        : data_(data)
        , size_(size)
    {
    }

    const_reference operator[](size_type i) const noexcept
    {
        return data_[i];
    }
    const_reference at(size_type i) const
    {
        check_if(i < size_, "index {} is out of range [0;{})", i, size_);
        return data_[i];
    }
    const_reference front() const
    {
        return data_[0];
    }
    const_reference back() const
    {
        return data_[size_ - 1];
    }
    size_type size() const
    {
        return size_;
    }
    const_pointer data() const noexcept
    {
        return data_;
    }

private:
    const_pointer data_;
    size_type size_;
};

#endif // GEOM_MODEL_SRC_UTIL_VECTOR_VIEW_HPP_