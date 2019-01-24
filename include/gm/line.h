#ifndef GEOM_MODEL_INCLUDE_GM_LINE_H_
#define GEOM_MODEL_INCLUDE_GM_LINE_H_

#include "abstract_curve.h"
#include "point.h"
#include "vec.h"

namespace gm {

class Line : public AbstractCurve {
public:
    Line() noexcept;
    Line(Vec dir, Point c) noexcept;

    const Point& c() const noexcept;
    const Vec& dir() const noexcept;

    Point f(double u) const noexcept override;
    Vec df(double u) const noexcept override;
    Vec df2(double u) const noexcept override;
    double project(const Point& x) const override;
    double approx_length(double begin, double end, size_t n) const override;

protected:
    std::ostream& print(std::ostream& os) const override;

private:
    Vec dir_;
    Point c_;
};

} // namespace gm

#endif // GEOM_MODEL_INCLUDE_GM_LINE_H_