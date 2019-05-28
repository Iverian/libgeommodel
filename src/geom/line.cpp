#include "C:\Users\trololo\Documents\Projects\libgeommodel\include\gm\vec.hpp"
#include <gm/line.hpp>
#include <util/debug.hpp>

#include <iostream>

#include <fmt/ostream.h>

namespace gm {

Line::Line() noexcept
    : dir_(1, 0, 0)
    , c_(0, 0, 0)
{
}

Line::Line(Vec dir, Point c) noexcept
    : dir_(std::move(dir))
    , c_(std::move(c))
{
}

const Point& Line::c() const noexcept
{
    return c_;
}

const Vec& Line::dir() const noexcept
{
    return dir_;
}

Point Line::f(double u) const noexcept
{
    return c_ + dir_ * u;
}

Vec Line::df(double u) const noexcept
{
    return dir_;
}

Vec Line::df2(double u) const noexcept
{
    return Vec();
}

double Line::project(const Point& p) const
{
    return dot(p - c_, dir_) / sqr(dir_);
}

double Line::approx_length(double begin, double end, size_t n) const
{
    return norm(dir_) * fabs(end - begin);
}

std::ostream& Line::print(std::ostream& os) const
{

    fmt::print(os, "{{ \"type\": \"line\", \"dir\": {0}, \"center\": {1} }}",
               dir_, c_);
    return os;
}

} // namespace gm