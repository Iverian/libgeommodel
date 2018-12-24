#ifndef GEOM_MODEL_INCLUDE_LINE_H_
#define GEOM_MODEL_INCLUDE_LINE_H_

#include "abstract_curve.h"

#include <memory>

class Line : public AbstractCurve {
public:
    ~Line() override;
    Line(Line&&) noexcept;
    Line& operator=(Line&&) noexcept;
    Line(const Line& other);
    Line& operator=(const Line& other);

    Line();
    Line(Vec dir, Point c);

    const Point& center() const;
    const Vec& dir() const;

    Point f(double u) const override;
    Vec df(double u) const override;
    Vec df2(double u) const override;
    double project(const Point& x) const override;
    double length(double begin, double end) const override;

protected:
    std::ostream& print(std::ostream& os) const override;

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl_;
};

#endif // GEOM_MODEL_INCLUDE_LINE_H_