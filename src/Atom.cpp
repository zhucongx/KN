#include "Atom.h"

#include <utility>

Atom::Atom(int id) : id_(id), mass_(0), type_("Vac") {}
Atom::Atom(int id, double mass, std::string type)
    : id_(id), mass_(mass), type_(std::move(type)) {}
Atom::Atom(int id, double mass, std::string type, double x, double y, double z)
    : id_(id), mass_(mass), type_(std::move(type)) {
  // Set both relative and absolute position, but will be corrected later
  relative_position_[kXDim] = x;
  relative_position_[kYDim] = y;
  relative_position_[kZDim] = z;

  absolute_position_[kXDim] = x;
  absolute_position_[kYDim] = y;
  absolute_position_[kZDim] = z;
}
Atom::~Atom() = default;
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

