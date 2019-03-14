#ifndef GEOM_MODEL_INCLUDE_GM_LINE_H_
#define GEOM_MODEL_INCLUDE_GM_LINE_H_

#include "abstract_curve.h"
#include "exports.h"
#include "point.h"
#include "vec.h"

namespace gm {

class GM_EXPORT Line : public AbstractCurve {
public:
    Line() noexcept;
    Line(Vec dir, Point c) noexcept;

    [[nodiscard]] const Point& c() const noexcept;
    [[nodiscard]] const Vec& dir() const noexcept;

    [[nodiscard]] Point f(double u) const noexcept override;
    [[nodiscard]] Vec df(double u) const noexcept override;
    [[nodiscard]] Vec df2(double u) const noexcept override;
    [[nodiscard]] double project(const Point& x) const override;
    [[nodiscard]] double approx_length(double begin, double end,
                                       size_t n) const override;

protected:
    std::ostream& print(std::ostream& os) const override;

private:
    Vec dir_;
    Point c_;
};

} // namespace gm

#endif // GEOM_MODEL_INCLUDE_GM_LINE_H_