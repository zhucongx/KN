#include "Atom.h"

namespace box {
Atom::Atom() = default;
Atom::Atom(Rank id, double mass, std::string type, double x, double y, double z)
    : id_(id), mass_(mass), type_(std::move(type)) {
  // Set both relative and absolute position, but will be corrected later
  relative_position_.x = x;
  relative_position_.y = y;
  relative_position_.z = z;

  absolute_position_.x = x;
  absolute_position_.y = y;
  absolute_position_.z = z;
}
Atom::Rank Atom::GetId() const {
  return id_;
}
void Atom::SetId(Rank id) {
  id_ = id;
}
double Atom::GetMass() const {
  return mass_;
}
void Atom::SetMass(double mass) {
  mass_ = mass;
}
const std::string &Atom::GetType() const {
  return type_;
}
void Atom::SetType(const std::string &type) {
  type_ = type;
}

Vector3<double> GetRelativeDistanceVector(const Atom &first, const Atom &second) {
  Vector3<double>
      relative_distance_vector = StarProduct(first.relative_position_,
                                             second.relative_position_);
  auto check_periodic = [](double &distance) {
    if (distance >= 0.5)
      distance -= 1;
    else if (distance < -0.5)
      distance += 1;
  };
  // periodic boundary conditions
  check_periodic(relative_distance_vector.x);
  check_periodic(relative_distance_vector.y);
  check_periodic(relative_distance_vector.z);
  return relative_distance_vector;
}
}// namespace box
