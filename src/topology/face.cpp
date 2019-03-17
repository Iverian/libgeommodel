#include <gm/face.hpp>

#include <util/util.hpp>

#include <fmt/ostream.h>

namespace gm {

struct Face::Impl {
    Impl(std::shared_ptr<AbstractSurface> surface, bool same_sense);
    Impl(std::shared_ptr<AbstractSurface> surface, bool same_sense,
         const FaceBound& outer_loop,
         const std::vector<FaceBound>& inner_loop);

    bool same_sense() const;
    const std::vector<FaceBound>& boundaries() const;
    const FaceBound& outer() const;
    const AbstractSurface& surface() const;

private:
    bool same_sense_;
    std::vector<FaceBound> boundaries_;
    std::shared_ptr<AbstractSurface> surface_;
};

Face::Face(std::shared_ptr<AbstractSurface> surface, bool same_sense)
    : pimpl_(std::make_shared<Face::Impl>(surface, same_sense))
{
}

Face::Face(std::shared_ptr<AbstractSurface> surface, bool same_sense,
           const FaceBound& outer_loop,
           const std::vector<FaceBound>& inner_loop)
    : pimpl_(std::make_shared<Face::Impl>(surface, same_sense, outer_loop,
                                          inner_loop))
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

const std::vector<FaceBound>& Face::boundaries() const
{
    return pimpl_->boundaries();
}

const AbstractSurface& Face::surface() const
{
    return pimpl_->surface();
}

std::ostream& operator<<(std::ostream& os, const Face& f)
{
    fmt::print(
        os, "{{ \"surface\": {0}, \"boundaries\": {1}, \"same_sense\": {2} }}",
        f.surface(), f.boundaries(), f.same_sense());
    return os;
}

Face::Impl::Impl(std::shared_ptr<AbstractSurface> surface, bool same_sense)
    : same_sense_(same_sense)
    , boundaries_()
    , surface_(move(surface))
{
}

Face::Impl::Impl(std::shared_ptr<AbstractSurface> surface, bool same_sense,
                 const FaceBound& outer_loop,
                 const std::vector<FaceBound>& inner_loop)
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

const std::vector<FaceBound>& Face::Impl::boundaries() const
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