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

Shell::Shell()
    : ax_()
    , faces_()
{
}

Shell::Shell(const Axis& ax, const vector<Face>& face)
    : ax_(ax)
    , faces_(face)
{
}

const vector<Face>& Shell::faces() const
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

void Shell::set_faces(const vector<Face>& faces)
{
    faces_ = faces;
}

void Shell::set_ax(Axis&& ax)
{
    ax_ = move(ax);
}

void Shell::set_faces(vector<Face>&& faces)
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