#include "Atom.h"

namespace kn
{
Atom::Atom() = default;
Atom::Atom(Rank id, double mass, std::string type, double x, double y, double z)
    : id_(id),
      mass_(mass),
      type_(std::move(type)),
      relative_position_{x, y, z},
      cartesian_position_{x, y, z}
{
}
Atom::Atom(Atom::Rank id,
           double mass,
           std::string type,
           Vector3 position)
    : id_(id),
      mass_(mass),
      type_(std::move(type)),
      relative_position_(position),
      cartesian_position_(position)
{
}
}// namespace kn
