#include <gm/oriented_edge.h>
#include <gm/shell.h>

#include <util/itertools.h>
#include <util/util.h>

#include <algorithm>
#include <optional>
#include <unordered_set>
#include <utility>

#include <fmt/ostream.h>

using namespace std;

namespace gm {

struct Shell::Impl {
    Impl();
    Impl(const Axis& ax, const vector<Face>& face);
    Impl(const Axis& ax, const vector<Face>& face, const vector<Edge>& edges);

    const vector<Face>& faces() const;
    const Axis& ax() const;
    const vector<Edge>& edges() const;

    void set_ax(const Axis& ax);
    void set_ax(Axis&& ax);
    void set_faces(const vector<Face>& faces);
    void set_faces(vector<Face>&& faces);

private:
    Axis ax_;
    vector<Face> faces_;
    mutable optional<vector<Edge>> edges_;
};

Shell::~Shell() = default;
Shell::Shell(Shell&&) noexcept = default;
Shell& Shell::operator=(Shell&&) noexcept = default;
Shell::Shell(const Shell& other)
    : pimpl_(make_unique<Impl>(*other.pimpl_))
{
}
Shell& Shell::operator=(const Shell& other)
{
    *pimpl_ = *other.pimpl_;
    return *this;
};

Shell::Shell()
    : pimpl_(make_unique<Impl>())
{
}

Shell::Shell(const Axis& ax, const vector<Face>& faces)
    : pimpl_(make_unique<Shell::Impl>(ax, faces))
{
}

Shell::Shell(const Axis& ax, const vector<Face>& faces,
             const vector<Edge>& edges)
    : pimpl_(make_unique<Shell::Impl>(ax, faces, edges))
{
}

const vector<Face>& Shell::faces() const
{
    return pimpl_->faces();
}

const vector<Edge>& Shell::edges() const
{
    return pimpl_->edges();
}

const Axis& Shell::ax() const
{
    return pimpl_->ax();
}

void Shell::set_ax(const Axis& ax)
{
    pimpl_->set_ax(ax);
}

void Shell::set_ax(Axis&& ax)
{
    pimpl_->set_ax(ax);
}

void Shell::set_faces(const vector<Face>& faces)
{
    pimpl_->set_faces(faces);
}

void Shell::set_faces(vector<Face>&& faces)
{
    pimpl_->set_faces(faces);
}

Shell::Impl::Impl()
    : ax_()
    , faces_()
{
}

Shell::Impl::Impl(const Axis& ax, const vector<Face>& face)
    : ax_(ax)
    , faces_(face)
{
}

Shell::Impl::Impl(const Axis& ax, const vector<Face>& face,
                  const vector<Edge>& edges)
    : ax_(ax)
    , faces_(face)
    , edges_(edges)
{
}

const vector<Face>& Shell::Impl::faces() const
{
    return faces_;
}

const Axis& Shell::Impl::ax() const
{
    return ax_;
}

const vector<Edge>& Shell::Impl::edges() const
{
    if (!edges_.has_value()) {
        unordered_set<Edge> unique;
        for (auto& f : faces_) {
            for (auto& i : f.boundaries()) {
                for (auto& j : i) {
                    unique.insert(j.edge());
                }
            }
        }
        edges_ = vector<Edge>(begin(unique), end(unique));
    }
    return edges_.value();
}

void Shell::Impl::set_ax(const Axis& ax)
{
    ax_ = ax;
}

void Shell::Impl::set_faces(const vector<Face>& faces)
{
    faces_ = faces;
}

void Shell::Impl::set_ax(Axis&& ax)
{
    ax_ = move(ax);
}

void Shell::Impl::set_faces(vector<Face>&& faces)
{
    faces_ = move(faces);
}

ostream& operator<<(ostream& os, const Shell& x)
{
    fmt::print(os, "{{ \"axis\": {0}, \"faces\": {1} }}", x.ax(),
               RangePrint(begin(x.faces()), end(x.faces())));
    return os;
}

} // namespace gm