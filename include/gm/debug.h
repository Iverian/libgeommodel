#ifndef GEOM_MODEL_INCLUDE_GM_DEBUG_H_
#define GEOM_MODEL_INCLUDE_GM_DEBUG_H_

namespace gm {

#ifdef NDEBUG

#define __GM_NOEXCEPT_RELEASE__ noexcept
static constexpr auto debug_flag = false;

#else // NDEBUG

#define __GM_NOEXCEPT_RELEASE__
static constexpr auto debug_flag = true;

#endif // NDEBUG

} // namespace gm

#endif // GEOM_MODEL_INCLUDE_GM_DEBUG_H_
