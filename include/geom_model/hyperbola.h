#ifndef GEOM_MODEL_INCLUDE_HYPERBOLA_H_
#define GEOM_MODEL_INCLUDE_HYPERBOLA_H_

#include "abstract_curve.h"
#include "axis.h"

#include <memory>

class Hyperbola : public AbstractCurve
{
  public:
    ~Hyperbola() override;
    Hyperbola(Hyperbola &&) noexcept;
    Hyperbola &operator=(Hyperbola &&) noexcept;
    Hyperbola(const Hyperbola &other);
    Hyperbola &operator=(const Hyperbola &other);

    Hyperbola();
    Hyperbola(double rx, double ry, Axis ax);

    Point f(double u) const override;
    Vec df(double u) const override;
    Vec df2(double u) const override;
    double project(const Point &p) const override;

  protected:
    std::ostream &print(std::ostream &os) const override;

  private:
    struct Impl;
    std::unique_ptr<Impl> pimpl_;
};

#endif // GEOM_MODEL_INCLUDE_HYPERBOLA_H_
