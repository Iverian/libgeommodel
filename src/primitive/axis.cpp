#include <geom_model/axis.h>

#include <fmt/ostream.h>

#include <util/math.h>

using namespace std;

struct Axis::Impl {
    Impl();
    Impl(const Vec& x, const Vec& y, const Vec& z, const Point& center);

    const Vec& at(size_t j) const;
    const Vec& operator[](size_t j) const;
    const Point& center() const;
    const Mat& get_transform() const;
    tuple<Point, Vec, Vec, Vec> get_view() const;
    double arg(const Point& p) const;

    Vec rotate_x(double angle, const Vec& v);
    Vec rotate_y(double angle, const Vec& v);
    Vec rotate_z(double angle, const Vec& v);

    Point rotate_x(double angle, const Point& p);
    Point rotate_y(double angle, const Point& p);
    Point rotate_z(double angle, const Point& p);

    Point global(const Point& p) const;
    Vec global(const Vec& v) const;

    bool operator==(const Impl& other) const;

private:
    Mat transform_;
    array<Vec, 3> basis_;
    Point c_;
};

Axis::~Axis() = default;
Axis::Axis(Axis&&) noexcept = default;
Axis& Axis::operator=(Axis&&) noexcept = default;
Axis::Axis(const Axis& other)
    : pimpl_(make_unique<Impl>(*other.pimpl_))
{
}
Axis& Axis::operator=(const Axis& other)
{
    *pimpl_ = *other.pimpl_;
    return *this;
};

Axis::Axis()
    : pimpl_(make_unique<Axis::Impl>())
{
}

Axis::Axis(const Vec& x, const Vec& y, const Vec& z, const Point& center)
    : pimpl_(make_unique<Axis::Impl>(x, y, z, center))
{
}

const Vec& Axis::at(size_t j) const
{
    return pimpl_->at(j);
}

const Vec& Axis::operator[](size_t j) const
{
    return (*pimpl_)[j];
}

const Vec& Axis::x() const
{
    return (*pimpl_)[0];
}

const Vec& Axis::y() const
{
    return (*pimpl_)[1];
}

const Vec& Axis::z() const
{
    return (*pimpl_)[2];
}

const Point& Axis::center() const
{
    return pimpl_->center();
}

const Mat& Axis::get_transform() const
{
    return pimpl_->get_transform();
}

tuple<Point, Vec, Vec, Vec> Axis::get_view() const
{
    return pimpl_->get_view();
}

Vec Axis::rotate_x(double angle, const Vec& v)
{
    return pimpl_->rotate_x(angle, v);
}

Vec Axis::rotate_y(double angle, const Vec& v)
{
    return pimpl_->rotate_y(angle, v);
}

Vec Axis::rotate_z(double angle, const Vec& v)
{
    return pimpl_->rotate_z(angle, v);
}

double Axis::arg(const Point& p) const
{
    return pimpl_->arg(p);
}

Point Axis::rotate_x(double angle, const Point& p)
{
    return pimpl_->rotate_x(angle, p);
}

Point Axis::rotate_y(double angle, const Point& p)
{
    return pimpl_->rotate_y(angle, p);
}

Point Axis::rotate_z(double angle, const Point& p)
{
    return pimpl_->rotate_z(angle, p);
}

Point Axis::global(const Point& p) const
{
    return pimpl_->global(p);
}

Vec Axis::global(const Vec& v) const
{
    return pimpl_->global(v);
}

bool operator==(const Axis& lhs, const Axis& rhs)
{
    return *lhs.pimpl_ == *rhs.pimpl_;
}

bool operator!=(const Axis& lhs, const Axis& rhs)
{
    return !(*lhs.pimpl_ == *rhs.pimpl_);
}

Point Axis::pglobal(double x, double y, double z) const
{
    return pimpl_->global(Point(x, y, z));
}

Vec Axis::vglobal(double x, double y, double z) const
{
    return pimpl_->global(Vec(x, y, z));
}

ostream& operator<<(ostream& os, const Axis& ax)
{
    fmt::print(os, "{{ \"x\": {0}, \"y\": {1}, \"z\": {2}, \"c\": {3} }}",
               ax[0], ax[1], ax[2], ax.center());
    return os;
}

Axis Axis::from_xy(const Vec& x, const Vec& ref_y, const Point& center)
{
    auto y = ref_y - dot(ref_y, x) * x;
    auto z = cross(x, y);
    return Axis(x, y, z, center);
}

Axis Axis::from_zx(const Vec& z, const Vec& ref_x, const Point& center)
{
    auto x = ref_x - dot(ref_x, z) * z;
    auto y = cross(z, x);
    return Axis(x, y, z, center);
}

Axis Axis::from_abc(const Point& a, const Point& b, const Point& c)
{
    return Axis::from_xy(Vec(a - b), Vec(c - b), b);
}

Axis::Impl::Impl()
    : Impl({1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0})
{
}

Axis::Impl::Impl(const Vec& x, const Vec& y, const Vec& z, const Point& center)
    : transform_()
    , basis_{x, y, z}
    , c_(center)
{
    ortogonalize(basis_);
    for (auto& i : basis_)
        i.normalize();
    transform_ = Mat({x[0], y[0], z[0], x[1], y[1], z[1], x[2], y[2], z[2]});
}

const Vec& Axis::Impl::at(size_t j) const
{
    return basis_.at(j);
}

const Vec& Axis::Impl::operator[](size_t j) const
{
    return basis_.at(j);
}

const Point& Axis::Impl::center() const
{
    return c_;
}

const Mat& Axis::Impl::get_transform() const
{
    return transform_;
}

tuple<Point, Vec, Vec, Vec> Axis::Impl::get_view() const
{
    return make_tuple<const Point&, const Vec&, const Vec&, const Vec&>(
        c_, basis_[0], basis_[1], basis_[2]);
}

double Axis::Impl::arg(const Point& p) const
{
    double px = dot(p - c_, basis_[0]), py = dot(p - c_, basis_[1]);
    return atan2(py, px);
}

Vec Axis::Impl::rotate_x(double angle, const Vec& v)
{
    return Mat::rotate(angle, basis_[0]) * v;
}

Vec Axis::Impl::rotate_y(double angle, const Vec& v)
{
    return Mat::rotate(angle, basis_[1]) * v;
}

Vec Axis::Impl::rotate_z(double angle, const Vec& v)
{
    return Mat::rotate(angle, basis_[2]) * v;
}

Point Axis::Impl::rotate_x(double angle, const Point& p)
{
    return Mat::rotate(angle, basis_[0]) * (p - c_) + c_;
}

Point Axis::Impl::rotate_y(double angle, const Point& p)
{
    return Mat::rotate(angle, basis_[1]) * (p - c_) + c_;
}

Point Axis::Impl::rotate_z(double angle, const Point& p)
{
    return Mat::rotate(angle, basis_[2]) * (p - c_) + c_;
}

Point Axis::Impl::global(const Point& p) const
{
    return dot(transform_, p) + c_;
}

Vec Axis::Impl::global(const Vec& v) const
{
    return dot(transform_, v);
}

bool Axis::Impl::operator==(const Axis::Impl& other) const
{
    return transform_ == other.transform_ && c_ == other.c_;
}
