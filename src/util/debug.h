#ifndef GEOM_MODEL_SRC_UTIL_DEBUG_H_
#define GEOM_MODEL_SRC_UTIL_DEBUG_H_

#include <iostream>
#include <stdexcept>

#include <fmt/ostream.h>

#ifdef NDEBUG

#define NOEXCEPTD noexcept
static constexpr auto debug_flag = false;

#else // NDEBUG

#define NOEXCEPTD
static constexpr auto debug_flag = true;

#endif // NDEBUG

#define DBG if (debug_flag)

#define cerrd                                                                 \
    if (!debug_flag) {                                                        \
    } else                                                                    \
        std::cerr

#define coutd                                                                 \
    if (!debug_flag) {                                                        \
    } else                                                                    \
        std::cout

#define throw_fmt(fmt_string, ...)                                            \
    do {                                                                      \
        auto what = fmt::format("({}:{}) " fmt_string, __FILE__, __LINE__,    \
                                ##__VA_ARGS__);                               \
        throw std::runtime_error(what);                                       \
    } while (0)

#define check_if(condition, fmt_string, ...)                                  \
    do {                                                                      \
        if (!(condition)) {                                                   \
            throw_fmt(fmt_string, ##__VA_ARGS__);                             \
        }                                                                     \
    } while (0)

#ifdef NDEBUG

#define debug_fmt(stream, fmt_string, ...)
#define check_ifd(condition, fmt_string, ...)

#else // NDEBUG

#define debug_fmt(stream, fmt_string, ...)                                    \
    do {                                                                      \
        fmt::print((stream), fmt_string, ##__VA_ARGS__);                      \
    } while (0)

#define check_ifd(condition, fmt_string, ...)                                 \
    check_if(condition, fmt_string, ##__VA_ARGS__)

#endif // NDEBUG

#endif // GEOM_MODEL_SRC_UTIL_DEBUG_H_
