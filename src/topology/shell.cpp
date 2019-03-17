#include <gm/oriented_edge.hpp>
#include <gm/shell.hpp>

#include <util/itertools.hpp>
#include <util/util.hpp>

#include <algorithm>
#include <optional>
#include <unordered_set>
#include <utility>

#include <fmt/ostream.h>

namespace gm {

Shell::Shell()
    : ax_()
    , faces_()
{
}

Shell::Shell(const Axis& ax, const std::vector<Face>& face)
    : ax_(ax)
    , faces_(face)
{
}

const std::vector<Face>& Shell::faces() const
{
    return faces_;
}

const Axis& Shell::ax() const
{
    return ax_;
}

void Shell::set_ax(const Axis& ax)
{
    ax_ = ax;
}

void Shell::set_faces(const std::vector<Face>& faces)
{
    faces_ = faces;
}

void Shell::set_ax(Axis&& ax)
{
    ax_ = std::move(ax);
}

void Shell::set_faces(std::vector<Face>&& faces)
{
    faces_ = std::move(faces);
}

std::ostream& operator<<(std::ostream& os, const Shell& x)
{
    fmt::print(os, "{{ \"axis\": {0}, \"faces\": {1} }}", x.ax(),
               RangePrint(std::begin(x.faces()), std::end(x.faces())));
    return os;
}

} // namespace gm