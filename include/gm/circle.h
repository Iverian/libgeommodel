#ifndef GEOM_MODEL_INCLUDE_CIRCLE_H_
#define GEOM_MODEL_INCLUDE_CIRCLE_H_

#include "axis.h"
#include "ellipse.h"

namespace gm {

class Circle : public Ellipse {
public:
    Circle() noexcept;
    Circle(double r, Axis ax) noexcept;

    double r() const noexcept;

    double approx_length(double begin, double end, size_t n) const override;

protected:
    std::ostream& print(std::ostream& os) const override;
};

} // namespace gm

#endif // GEOM_MODEL_INCLUDE_CIRCLE_H_