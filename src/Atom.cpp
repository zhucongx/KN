#include "Atom.h"

namespace box
{
Atom::Atom() = default;
Atom::Atom(Rank id, double mass, std::string type, double x, double y, double z)
    : id_(id), mass_(mass), type_(std::move(type))
{
  // Set both relative and absolute position, but will be corrected later
  relative_position_[kXDimension] = x;
  relative_position_[kYDimension] = y;
  relative_position_[kZDimension] = z;

  cartesian_position_[kXDimension] = x;
  cartesian_position_[kYDimension] = y;
  cartesian_position_[kZDimension] = z;
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

Atom::Rank Atom::GetId() const
{
  return id_;
}
void Atom::SetId(Rank id)
{
  id_ = id;
}
double Atom::GetMass() const
{
  return mass_;
}
void Atom::SetMass(double mass)
{
  mass_ = mass;
}
const std::string &Atom::GetType() const
{
  return type_;
}
void Atom::SetType(const std::string &type)
{
  type_ = type;
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

}// namespace box
