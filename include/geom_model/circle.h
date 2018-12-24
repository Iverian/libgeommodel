#ifndef GEOM_MODEL_INCLUDE_CIRCLE_H_
#define GEOM_MODEL_INCLUDE_CIRCLE_H_

#include "axis.h"
#include "ellipse.h"

#include <memory>

class Circle : public Ellipse {
public:
    ~Circle() override;
    Circle(Circle&&) noexcept;
    Circle& operator=(Circle&&) noexcept;
    Circle(const Circle& other);
    Circle& operator=(const Circle& other);

    Circle();
    Circle(double r, Axis ax);

    double r() const;

    double length(double begin, double end) const override;

protected:
    std::ostream& print(std::ostream& os) const override;

private:
    struct Impl;
    std::unique_ptr<Circle::Impl> pimpl_;
};

#endif // GEOM_MODEL_INCLUDE_CIRCLE_H_