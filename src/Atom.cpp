#include "Atom.h"
namespace kn {
Atom::Atom() = default;

Atom::Atom(int id, double mass, std::string type, double x, double y, double z) :
    cartesian_position_{x, y, z},
    relative_position_{x, y, z},
    id_(id),
    mass_(mass),
    type_(std::move(type)) {
  first_nearest_neighbor_list_.reserve(Al_const::kNumFirstNearestNeighbors);
  near_neighbor_list_.reserve(Al_const::kNumNearNeighbors);
}

Atom::Atom(int id, double mass, std::string type, Vector3 position) :
    cartesian_position_(position),
    relative_position_(position),
    id_(id),
    mass_(mass),
    type_(std::move(type)) {
  first_nearest_neighbor_list_.reserve(Al_const::kNumFirstNearestNeighbors);
  near_neighbor_list_.reserve(Al_const::kNumNearNeighbors);
}
const Vector3 &Atom::GetCartesianPosition() const {
  return cartesian_position_;
}
void Atom::SetCartesianPosition(const Vector3 &cartesian_position) {
  cartesian_position_ = cartesian_position;
}
const Vector3 &Atom::GetRelativePosition() const {
  return relative_position_;
}

void Atom::SetRelativePosition(const Vector3 &relative_position) {
  relative_position_ = relative_position;
}
const std::vector<int> &Atom::GetNearNeighborList() const {
  return near_neighbor_list_;
}
const std::vector<int> &Atom::GetFirstNearestNeighborList() const {
  return first_nearest_neighbor_list_;
}

int Atom::GetId() const {
  return id_;
}
double Atom::GetMass() const {
  return mass_;
}
const std::string &Atom::GetType() const {
  return type_;
}
void Atom::SetType(const std::string &type) {
  type_ = type;
}
void Atom::AppendNearNeighborList(int index) {
  near_neighbor_list_.emplace_back(index);
}
void Atom::AppendFirstNearestNeighborList(int index) {
  first_nearest_neighbor_list_.emplace_back(index);
}
Vector3 GetRelativeDistanceVector(const Atom &first, const Atom &second) {
  Vector3 relative_distance_vector = first.relative_position_ - second.relative_position_;
  auto check_periodic = [](double &distance) {
    if (distance >= 0.5)
      distance -= 1;
    else if (distance < -0.5)
      distance += 1;
  };
  // periodic boundary conditions
  for (const auto kDim : All_Dimensions) {
    check_periodic(relative_distance_vector[kDim]);
  }
  return relative_distance_vector;
}


} // namespace kn