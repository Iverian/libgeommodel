#ifndef GEOM_MODEL_INCLUDE_GM_CIRCLE_H_
#define GEOM_MODEL_INCLUDE_GM_CIRCLE_H_

#include "axis.h"
#include "ellipse.h"
#include "exports.h"

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

#endif // GEOM_MODEL_INCLUDE_GM_CIRCLE_H_