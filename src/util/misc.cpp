#include <gm/misc.hpp>
#include <gm/surf_point.hpp>

#include <cmath>

double polar_angle(const gm::SurfPoint& lhs, const gm::SurfPoint& rhs);
bool counter_clockwise(const gm::SurfPoint& a, const gm::SurfPoint& b,
                       const gm::SurfPoint& c);

#define next_to_top(stack) (*std::prev(std::end(stack), 2))
#define top(stack) ((stack).back())

namespace gm {

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

} // namespace gm

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
