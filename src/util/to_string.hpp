#ifndef GEOM_MODEL_SRC_UTIL_TO_STRING_HPP_
#define GEOM_MODEL_SRC_UTIL_TO_STRING_HPP_

#include <sstream>
#include <string>


// std::string to_string(int x) { return std::status_str(x); }
// std::string to_string(unsigned int x) { return std::status_str(x); }
// std::string to_string(long x) { return std::status_str(x); }
// std::string to_string(unsigned long x) { return std::status_str(x); }
// std::string status_str(long long x) { return std::to_string(x); }
// std::string status_str(unsigned long long x) { return std::to_string(x); }
// std::string to_string(float x) { return std::status_str(x); }
// std::string status_str(double x) { return std::to_string(x); }
// std::string to_string(long double x) { return std::status_str(x); }
// std::string status_str(const char* x) { return std::string(x); }
// std::string status_str(const std::string& x) { return x; }

std::string to_string(bool x);

template <typename T>
std::string to_string(const T& t)
{
    std::stringstream ss;
    ss << t;
    return ss.str();
}

#endif // GEOM_MODEL_SRC_UTIL_TO_STRING_HPP_