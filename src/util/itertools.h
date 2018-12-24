#ifndef GEOM_MODEL_SRC_UTIL_ITERTOOLS_H_
#define GEOM_MODEL_SRC_UTIL_ITERTOOLS_H_

#include <functional>
#include <iostream>
#include <iterator>
#include <optional>

template <class ForwardIt>
struct RangePrint {
    RangePrint(ForwardIt first, ForwardIt last)
        : first_(first)
        , last_(last)
    {
    }

    ForwardIt first() const
    {
        return first_;
    }

    ForwardIt last() const
    {
        return last_;
    }

private:
    ForwardIt first_;
    ForwardIt last_;
};

template <class ForwardIt>
std::ostream& operator<<(std::ostream& os, const RangePrint<ForwardIt>& obj)
{
    os << "[";
    for (auto i = obj.first(); i != obj.last(); ++i) {
        os << (*i) << ((std::next(i) != obj.last()) ? ", " : "]");
    }
    return os;
}

#endif // GEOM_MODEL_SRC_UTIL_ITERTOOLS_H_