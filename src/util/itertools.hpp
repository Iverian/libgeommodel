#ifndef GEOM_MODEL_SRC_UTIL_ITERTOOLS_HPP_
#define GEOM_MODEL_SRC_UTIL_ITERTOOLS_HPP_

#include <functional>
#include <iostream>
#include <iterator>
#include <optional>

template <class ForwardIt>
struct RangePrint {
    RangePrint(ForwardIt pfirst, ForwardIt plast)
        : first(std::move(pfirst))
        , last(std::move(plast))
    {
    }

    template <class T>
    RangePrint(const T& container)
        : first(std::begin(container))
        , last(std::end(container))
    {
    }

    ForwardIt first;
    ForwardIt last;
};

template <class ForwardIt>
std::ostream& operator<<(std::ostream& os, const RangePrint<ForwardIt>& obj)
{
    os << "[";
    for (auto i = obj.first; i != obj.last; ++i) {
        os << (*i) << ((std::next(i) != obj.last) ? ", " : "]");
    }
    return os;
}

#endif // GEOM_MODEL_SRC_UTIL_ITERTOOLS_HPP_