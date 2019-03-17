#ifndef GEOM_MODEL_INCLUDE_GM_SHELL_HPP_
#define GEOM_MODEL_INCLUDE_GM_SHELL_HPP_

#include "axis.h"
#include "exports.h"
#include "face.h"
#include "oriented_edge.h"

namespace gm {

class GM_EXPORT Shell {
public:
    Shell();
    Shell(const Axis& ax, const std::vector<Face>& faces);

    [[nodiscard]] const std::vector<Face>& faces() const;
    [[nodiscard]] const Axis& ax() const;

    void set_ax(const Axis& ax);
    void set_ax(Axis&& ax);
    void set_faces(const std::vector<Face>& faces);
    void set_faces(std::vector<Face>&& faces);

private:
    Axis ax_;
    std::vector<Face> faces_;
};

GM_EXPORT std::ostream& operator<<(std::ostream& os, const Shell& x);

} // namespace gm

#endif // GEOM_MODEL_INCLUDE_GM_SHELL_HPP_
