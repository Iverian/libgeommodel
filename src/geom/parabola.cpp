#define _USE_MATH_DEFINES

#include <gm/compare.hpp>
#include <gm/parabola.hpp>

#include <util/math.hpp>

#include <fmt/ostream.hpp>

using namespace std;

namespace gm {

Parabola::Parabola() noexcept
    : f_(1)
    , ax_()
{
}

Parabola::Parabola(double f, Axis ax) noexcept
    : f_(fabs(f))
    , ax_(move(ax))
{
}

Point Parabola::f(double u) const noexcept
{
    return ax_.pglobal(f_ * ::sqr(u), 2 * u * f_, 0);
}

Vec Parabola::df(double u) const noexcept
{
    return ax_.vglobal(2 * f_ * u, 2 * f_, 0);
}

Vec Parabola::df2(double u) const noexcept
{
    return ax_.vglobal(2 * f_, 0, 0);
}

ostream& Parabola::print(ostream& os) const
{
    fmt::print(os, "{{ \"type\": \"parabola\", \"f\": {0}, \"axis\": {1} }}",
               f_, ax_);
    return os;
}

// проекция: решение уравнения u^3 + P u + Q == 0
double Parabola::project(const Point& p) const
{
    auto [c, x, y, z] = ax_.get_view();
    auto w = p - c;

    vector<double> roots;
    roots.reserve(3);

    auto P = 2 - dot(w, x) / (2 * f_);
    auto Q = -dot(w, y) / (2 * f_);
    auto D = ::pow(P / 3, 3) - ::sqr(Q / 2);

    if (cmp::zero(D)) {
        auto x = ::pow(-Q / 2, 1. / 3);
        if (cmp::zero(x)) {
            roots.emplace_back(0);
        } else {
            roots.emplace_back(x);
            roots.emplace_back(2 * x);
        }
    } else if (D > 0) {
        auto y = ::sqrt(D);
        auto x = ::pow(-Q / 2 + y, 1. / 3) + ::pow(-Q / 2 - y, 1. / 3);
        roots.emplace_back(x);
    } else {
        auto A = ::sqrt(-P / 3);
        auto B = ::acos(3 * Q * A / (2 * P)) / 3;
        auto C = 2 * M_PI / 3;

        for (auto i = 0; i < 3; ++i) {
            roots.emplace_back(2 * A * ::cos(B - C * i));
        }
    }

    return *min_element(begin(roots), end(roots),
                        [this, &p](auto& lhs, auto& rhs) {
                            return dist(f(lhs), p) < dist(f(rhs), p);
                        });
}

inline double F(double x)
{
    return x * ::sqrt(1 + ::sqr(x)) + ::asinh(x);
}

double Parabola::approx_length(double begin, double end, size_t n) const
{
    if (begin >= end) {
        swap(begin, end);
    }
    return 2 * f_ * (F(end) - F(begin));
}

} // namespace gm