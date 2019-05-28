#include "util.hpp"

#include <cmms/cyclic_iterator.hpp>
#include <gm/compare.hpp>
#include <util/vector_view.hpp>

std::optional<double>
zero_intersect(const std::pair<gm::SurfPoint, gm::SurfPoint>& line)
{
    auto& [a, b] = line;

    if (gm::cmp::zero(a.v, gm::Tolerance::MAX)) {
        return a.u;
    }
    if (gm::cmp::zero(b.v, gm::Tolerance::MAX)) {
        return b.u;
    }
    if (std::signbit(a.v) != std::signbit(b.v)) {
        return (a.u * b.v - b.u * a.v) / (b.v - a.v);
    }

    return std::nullopt;
}

std::vector<double>
single_eliminate(const std::vector<gm::SurfPoint>& convex_hull, double pfront,
                 double pback, gm::Tolerance tol) noexcept
{
    static constexpr auto max_roots = size_t(2);

    std::vector<double> roots;
    roots.reserve(max_roots);

    auto it = cmms::make_cycle(convex_hull);
    do {
        if (auto u = zero_intersect({*it, *std::next(it)}); u) {
            auto v = u.value();
            if (!gm::cmp::near(v, pfront, tol) && !gm::cmp::near(v, pback, tol)
                && (roots.empty() || !gm::cmp::near(v, roots.back()))) {
                roots.emplace_back(v);
            }
        }
    } while (++it, roots.size() != max_roots && it.iter() != it.first());
    std::sort(std::begin(roots), std::end(roots));

    return roots;
}

size_t find_span(double t, size_t order,
                 const VectorView<double>& knots) noexcept
{
    auto result = size_t(-1);
    auto last = knots.size() - order;

    if (gm::cmp::near(t, knots.back())) {
        result = last - 1;
    } else {
        for (auto i = order - 1; i <= last; ++i) {
            if (t < knots[i]) {
                result = i - 1;
                break;
            }
        }
    }

    return result;
}

double pget(const gm::WPoint<double, 1>& cp) noexcept
{
    return cp.p()[0];
}

double bord_check(double r, double a, double b)
{
    if (r < a) {
        r = a;
    }
    if (r > b) {
        r = b;
    }
    return r;
}

double bord_bounce(double r, double a, double b)
{
    while (r < a || r > b) {
        if (r < a) {
            r = 2 * a - r;
        }
        if (r > b) {
            r = 2 * b - r;
        }
    }
    return r;
}