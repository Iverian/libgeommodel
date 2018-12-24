#ifndef GEOM_MODEL_INCLUDE_ELLIPSE_H_
#define GEOM_MODEL_INCLUDE_ELLIPSE_H_

#include "abstract_curve.h"
#include "axis.h"

#include <memory>

class Ellipse : public AbstractCurve {
public:
    ~Ellipse() override;
    Ellipse(Ellipse&&) noexcept;
    Ellipse& operator=(Ellipse&&) noexcept;
    Ellipse(const Ellipse& other);
    Ellipse& operator=(const Ellipse& other);

    Ellipse();
    Ellipse(double rx, double ry, Axis ax);

    double rx() const;
    double ry() const;
    const Axis& ax() const;

    Point f(double u) const override;
    Vec df(double u) const override;
    Vec df2(double u) const override;
    double project(const Point& p) const override;
    std::optional<double> project_greater(const Point& p,
                                          double min) const override;

protected:
    std::ostream& print(std::ostream& os) const override;

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl_;
};

#endif // GEOM_MODEL_INCLUDE_ELLIPSE_H_
