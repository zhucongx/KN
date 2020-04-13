#include "Atom.h"

#include <utility>
namespace box {
Atom::Atom(int id, double mass, std::string type, double x, double y, double z)
    : id_(id), mass_(mass), type_(std::move(type)) {
  // Set both relative and absolute position, but will be corrected later
  relative_position_.x = x;
  relative_position_.y = y;
  relative_position_.z = z;

  absolute_position_.x = x;
  absolute_position_.y = y;
  absolute_position_.z = z;
}
int Atom::GetId() const {
  return id_;
}
void Atom::SetId(int id) {
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

void Atom::SetAbsolutePosition(
    const Double3 &absolute_position) {
  absolute_position_ = absolute_position;
}
void Atom::SetRelativePosition(
    const Double3 &relative_position) {
  relative_position_ = relative_position;
}
const Double3 &Atom::GetAbsolutePosition() const {
  return absolute_position_;
}
const Double3 &Atom::GetRelativePosition() const {
  return relative_position_;
}

Double3 GetRelativeDistanceVector(const Atom& first, const Atom& second){
  auto atom1_relative_position = first.GetRelativePosition();
  auto atom2_relative_position = second.GetRelativePosition();
  Double3 relative_distance_vector =
      {atom2_relative_position.x - atom1_relative_position.x,
       atom2_relative_position.y - atom1_relative_position.y,
       atom2_relative_position.z - atom1_relative_position.z};

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
