#ifndef GEOM_MODEL_INCLUDE_PARABOLA_H_
#define GEOM_MODEL_INCLUDE_PARABOLA_H_

#include "abstract_curve.h"
#include "axis.h"

#include <memory>

class Parabola : public AbstractCurve {
public:
    ~Parabola() override;
    Parabola(Parabola&&) noexcept;
    Parabola& operator=(Parabola&&) noexcept;
    Parabola(const Parabola& other);
    Parabola& operator=(const Parabola& other);

    Parabola();
    Parabola(double f, Axis ax);

    Point f(double u) const override;
    Vec df(double u) const override;
    Vec df2(double u) const override;
    double project(const Point& p) const override;

protected:
    std::ostream& print(std::ostream& os) const override;

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl_;
};

#endif // GEOM_MODEL_INCLUDE_PARABOLA_H_
