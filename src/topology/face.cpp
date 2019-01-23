#include <gm/face.h>

#include <util/util.h>

#include <fmt/ostream.h>

using namespace std;

namespace gm {

struct Face::Impl {
    Impl(shared_ptr<AbstractSurface> surface, bool same_sense);
    Impl(shared_ptr<AbstractSurface> surface, bool same_sense,
         const FaceBound& outer_loop, const vector<FaceBound>& inner_loop);

    bool same_sense() const;
    const vector<FaceBound>& boundaries() const;
    const FaceBound& outer() const;
    const AbstractSurface& surface() const;

private:
    bool same_sense_;
    vector<FaceBound> boundaries_;
    shared_ptr<AbstractSurface> surface_;
};

Face::Face(shared_ptr<AbstractSurface> surface, bool same_sense)
    : pimpl_(make_shared<Face::Impl>(surface, same_sense))
{
}

Face::Face(shared_ptr<AbstractSurface> surface, bool same_sense,
           const FaceBound& outer_loop, const vector<FaceBound>& inner_loop)
    : pimpl_(
          make_shared<Face::Impl>(surface, same_sense, outer_loop, inner_loop))
{
}

bool Face::same_sense() const
{
    return pimpl_->same_sense();
}

const FaceBound& Face::outer() const
{
    return pimpl_->outer();
}

const vector<FaceBound>& Face::boundaries() const
{
    return pimpl_->boundaries();
}

const AbstractSurface& Face::surface() const
{
    return pimpl_->surface();
}

ostream& operator<<(ostream& os, const Face& f)
{
    fmt::print(
        os, "{{ \"surface\": {0}, \"boundaries\": {1}, \"same_sense\": {2} }}",
        f.surface(), f.boundaries(), f.same_sense());
    return os;
}

Face::Impl::Impl(shared_ptr<AbstractSurface> surface, bool same_sense)
    : same_sense_(same_sense)
    , boundaries_()
    , surface_(move(surface))
{
}

Face::Impl::Impl(shared_ptr<AbstractSurface> surface, bool same_sense,
                 const FaceBound& outer_loop,
                 const vector<FaceBound>& inner_loop)
    : same_sense_(same_sense)
    , boundaries_()
    , surface_(move(surface))
{
    boundaries_.reserve(1 + inner_loop.size());
    boundaries_.emplace_back(outer_loop);
    for (auto& i : inner_loop) {
        boundaries_.emplace_back(i);
    }
}

const vector<FaceBound>& Face::Impl::boundaries() const
{
    return boundaries_;
}

const FaceBound& Face::Impl::outer() const
{
    return boundaries_.front();
}

bool Face::Impl::same_sense() const
{
    return same_sense_;
}

const AbstractSurface& Face::Impl::surface() const
{
    return *surface_;
}

} // namespace gm