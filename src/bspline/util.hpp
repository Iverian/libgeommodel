#ifndef GEOM_MODEL_SRC_BSPLINE_UTIL_HPP_
#define GEOM_MODEL_SRC_BSPLINE_UTIL_HPP_

#include <gm/surf_point.hpp>
#include <util/vector_view.hpp>

#include <optional>
#include <utility>
#include <vector>

std::vector<gm::SurfPoint> graham_scan(std::vector<gm::SurfPoint>& points);
double polar_angle(const gm::SurfPoint& lhs, const gm::SurfPoint& rhs);
bool counter_clockwise(const gm::SurfPoint& a, const gm::SurfPoint& b,
                       const gm::SurfPoint& c);
std::optional<double>
zero_intersect(const std::pair<gm::SurfPoint, gm::SurfPoint>& line);
std::vector<double>
single_eliminate(const std::vector<gm::SurfPoint>& convex_hull, double pfront,
                 double pback) noexcept;
size_t find_span(double t, size_t order,
                 const VectorView<double>& knots) noexcept;

#endif // GEOM_MODEL_SRC_BSPLINE_UTIL_HPP_