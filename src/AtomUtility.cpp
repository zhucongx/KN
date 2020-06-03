#include "AtomUtility.h"
namespace kn {
Vector3 AtomUtility::GetRelativeDistanceVector(const Atom &first, const Atom &second) {
  Vector3 relative_distance_vector = first.relative_position_ - second.relative_position_;
  auto check_periodic = [](double &distance) {
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