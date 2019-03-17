#ifndef GEOM_MODEL_SRC_UTIL_MATH_HPP_
#define GEOM_MODEL_SRC_UTIL_MATH_HPP_

#include <gm/surf_point.hpp>

#include "debug.h"
#include "std_defines.h"

#include <array>
#include <cmath>
#include <vector>

#include <type_traits>

double binom(size_t n, size_t k) __GM_NOEXCEPT_RELEASE__;

double pad(double t, double f, double b, double pad = 1e-9);
std::array<double, 2> atan2v(double y, double x);

template <class T>
T sqr(T x)
{
    return T(x * x);
}

template <class T, class Function>
T diff(Function f, double t)
{
    double h = 1e-4;
    return (-3 * f(t) + 4 * f(t + h) - f(t + 2 * h)) / (2 * h);
}

template <class T, class Function>
T diff2(Function f, double t)
{
    double h = 1e-2;
    return (-f(t + 2 * h) + f(t + h) * 16 - f(t) * 30 + f(t - h) * 16
            - f(t - 2 * h))
        / (12 * sqr(h));
}

template <class T, class Function>
T diff11(Function f, gm::SurfPoint t)
{
    double h = 1e-4;
    return (-3
                * (-3 * f({t.u, t.v}) + 4 * f({t.u, t.v + h})
                   - f({t.u, t.v + 2 * h}))
            + 4
                * (-3 * f({t.u + h, t.v}) + 4 * f({t.u + h, t.v + h})
                   - f({t.u + h, t.v + 2 * h}))
            - (-3 * f({t.u + 2 * h, t.v}) + 4 * f({t.u + 2 * h, t.v + h})
               - f({t.u + 2 * h, t.v + 2 * h})))
        / (4 * sqr(h));
}

template <class BidirIt,
          class = std::enable_if_t<
              std::is_floating_point_v<typename BidirIt::value_type>>>
typename BidirIt::value_type trapz(BidirIt first, BidirIt last,
                                   typename BidirIt::value_type step)
{
    typename BidirIt::value_type result {};

    typename BidirIt::value_type prev, cur = *first;
    for (auto i = std::next(first); i != last; ++i) {
        prev = cur;
        cur = *i;
        result += step * (prev + cur) / 2;
    }

    return result;
}

template <class Integer>
Integer floor(double x)
{
    return Integer(std::floor(x));
}

template <class Integer>
Integer ceil(double x)
{
    return Integer(std::ceil(x));
}

#endif // GEOM_MODEL_SRC_UTIL_MATH_HPP_
