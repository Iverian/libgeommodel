#ifndef GEOM_MODEL_INCLUDE_PLANE_H_
#define GEOM_MODEL_INCLUDE_PLANE_H_

#include "abstract_surface.h"
#include "axis.h"

#include <memory>

class Plane : public AbstractSurface {
public:
    ~Plane() override;
    Plane(Plane&&) noexcept;
    Plane& operator=(Plane&&) noexcept;
    Plane(const Plane& other);
    Plane& operator=(const Plane& other);

    Plane();
    explicit Plane(Axis ax);

    const Axis& ax() const;

    Point f(const ParametricPoint& p) const override;
    Vec dfu(const ParametricPoint& p) const override;
    Vec dfv(const ParametricPoint& p) const override;
    Vec dfuu(const ParametricPoint& p) const override;
    Vec dfvv(const ParametricPoint& p) const override;
    Vec dfuv(const ParametricPoint& p) const override;
    std::ostream& print(std::ostream& os) const override;
    ParametricPoint project(const Point& p) const override;

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl_;
};

#endif // GEOM_MODEL_INCLUDE_PLANE_H_
