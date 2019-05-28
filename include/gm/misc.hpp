#ifndef GEOM_MODEL_INCLUDE_GM_MISC_HPP_
#define GEOM_MODEL_INCLUDE_GM_MISC_HPP_

#include "surf_point.hpp"

#include <vector>

namespace gm {
std::vector<gm::SurfPoint> graham_scan(std::vector<gm::SurfPoint>& points);
} // namespace gm

#endif // GEOM_MODEL_INCLUDE_GM_MISC_HPP_