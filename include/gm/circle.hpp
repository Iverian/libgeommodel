#ifndef GEOM_MODEL_INCLUDE_GM_CIRCLE_HPP_
#define GEOM_MODEL_INCLUDE_GM_CIRCLE_HPP_

#include "axis.hpp"
#include "ellipse.hpp"
#include "exports.hpp"

namespace gm {

class GM_EXPORT Circle : public Ellipse {
public:
    Circle() noexcept;
    Circle(double r, Axis ax) noexcept;

    [[nodiscard]] double r() const noexcept;

    [[nodiscard]] double approx_length(double begin, double end,
                                       size_t n) const override;

protected:
    std::ostream& print(std::ostream& os) const override;
};

} // namespace gm

#endif // GEOM_MODEL_INCLUDE_GM_CIRCLE_HPP_