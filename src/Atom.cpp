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


Vector3 GetRelativeDistanceVector(const Atom &first, const Atom &second)
{
  Vector3 relative_distance_vector = first.relative_position_ - second.relative_position_;
  auto check_periodic = [](double &distance)
  {
    if (distance >= 0.5)
      distance -= 1;
    else if (distance < -0.5)
      distance += 1;
  };
  // periodic boundary conditions
  check_periodic(relative_distance_vector[kXDimension]);
  check_periodic(relative_distance_vector[kYDimension]);
  check_periodic(relative_distance_vector[kZDimension]);
  return relative_distance_vector;
}

}// namespace kn
