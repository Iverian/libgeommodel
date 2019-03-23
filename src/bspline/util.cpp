#include "util.hpp"

#include <gm/compare.hpp>
#include <util/cyclic_iterator.hpp>
#include <util/vector_view.hpp>

#define next_to_top(stack) (*std::prev(std::end(stack), 2))
#define top(stack) ((stack).back())

std::vector<gm::SurfPoint> graham_scan(std::vector<gm::SurfPoint>& points)
{
    std::vector<gm::SurfPoint> result;
    result.reserve(points.size());

    {
        auto min = std::min_element(
            std::begin(points), std::end(points),
            [](auto& lhs, auto& rhs) { return lhs.v < rhs.v; });
        std::swap(*std::begin(points), *min);
    }
    std::sort(std::next(std::begin(points)), std::end(points),
              [& p = points[0]](auto& lhs, auto& rhs) {
                  return polar_angle(lhs, p) < polar_angle(rhs, p);
              });
    result.push_back(points[0]);
    result.push_back(points[1]);

    for (size_t i = 2; i < points.size(); ++i) {
        auto& q = points[i];
        while (result.size() > 1
               && !counter_clockwise(next_to_top(result), top(result), q)) {
            result.pop_back();
        }
        result.push_back(q);
    }

    result.shrink_to_fit();
    return result;
}

#undef next_to_top
#undef top

double polar_angle(const gm::SurfPoint& lhs, const gm::SurfPoint& rhs)
{
    auto p = rhs - lhs;
    return std::atan2(p.v, p.u);
}

bool counter_clockwise(const gm::SurfPoint& a, const gm::SurfPoint& b,
                       const gm::SurfPoint& c)
{
    auto value = (b.u - a.u) * (c.v - a.v) - (b.v - a.v) * (c.u - a.u);
    return value > 0;
}

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
                 double pback) noexcept
{
    static constexpr auto max_roots = size_t(2);

    std::vector<double> roots;
    roots.reserve(max_roots);
    auto it = CyclicIterator(std::begin(convex_hull), std::end(convex_hull));
    do {
        if (auto u = zero_intersect({*it, *std::next(it)}); u) {
            auto v = u.value();
            if (!gm::cmp::near(v, pfront) && !gm::cmp::near(v, pback)
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