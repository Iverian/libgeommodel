#ifndef GEOM_MODEL_INCLUDE_SHELL_H_
#define GEOM_MODEL_INCLUDE_SHELL_H_

#include "axis.h"
#include "face.h"
#include "oriented_edge.h"

namespace gm {

struct Shell {
    ~Shell();
    Shell(Shell&&) noexcept;
    Shell& operator=(Shell&&) noexcept;
    Shell(const Shell& other);
    Shell& operator=(const Shell& other);

    Shell();
    Shell(const Axis& ax, const std::vector<Face>& faces);
    Shell(const Axis& ax, const std::vector<Face>& faces,
          const std::vector<Edge>& edges);

    const std::vector<Face>& faces() const;
    const Axis& ax() const;
    const std::vector<Edge>& edges() const;

    void set_ax(const Axis& ax);
    void set_ax(Axis&& ax);
    void set_faces(const std::vector<Face>& faces);
    void set_faces(std::vector<Face>&& faces);

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl_;
};

using BoundaryRep = std::vector<Shell>;

std::ostream& operator<<(std::ostream& os, const Shell& x);

} // namespace gm

#endif // GEOM_MODEL_INCLUDE_SHELL_H_
