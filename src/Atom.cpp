#include "Atom.h"

namespace kn {
Atom::Atom() = default;

Atom::Atom(Rank id, double mass, std::string type, double x, double y, double z) :
  cartesian_position_{x, y, z},
  relative_position_{x, y, z},
  id_(id),
  mass_(mass),
  type_(std::move(type)) {
}

Atom::Atom(Rank id,
           double mass,
           std::string type,
           Vector3 position) :
  cartesian_position_(position),
  relative_position_(position),
  id_(id),
  mass_(mass),
  type_(std::move(type)) {
}
} // namespace kn
