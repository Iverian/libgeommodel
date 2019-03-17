#ifndef GEOM_MODEL_SRC_BSPLINE_UTIL_HPP_
#define GEOM_MODEL_SRC_BSPLINE_UTIL_HPP_

#include <gm/surf_point.hpp>

#include <optional>
#include <utility>
#include <vector>

namespace gm {

std::vector<gm::SurfPoint> graham_scan(std::vector<gm::SurfPoint>& points);
double polar_angle(const gm::SurfPoint& lhs, const gm::SurfPoint& rhs);
bool counter_clockwise(const gm::SurfPoint& a, const gm::SurfPoint& b,
                       const gm::SurfPoint& c);
std::optional<double>
zero_intersect(const std::pair<gm::SurfPoint, gm::SurfPoint>& line);

} // namespace gm

#endif // GEOM_MODEL_SRC_BSPLINE_UTIL_HPP_