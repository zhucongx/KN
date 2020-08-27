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
  second_nearest_neighbor_list_.reserve(Al_const::kNumSecondNearestNeighbors);
  third_nearest_neighbor_list_.reserve(Al_const::kNumThirdNearestNeighbors);
}
Atom::Atom(int id, double mass, std::string type, Vector3 position) :
    cartesian_position_(position),
    relative_position_(position),
    id_(id),
    mass_(mass),
    type_(std::move(type)) {
  first_nearest_neighbor_list_.reserve(Al_const::kNumFirstNearestNeighbors);
  second_nearest_neighbor_list_.reserve(Al_const::kNumSecondNearestNeighbors);
  third_nearest_neighbor_list_.reserve(Al_const::kNumThirdNearestNeighbors);
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
void Atom::AppendFirstNearestNeighborList(int index) {
  first_nearest_neighbor_list_.emplace_back(index);
}
void Atom::AppendSecondNearestNeighborList(int index) {
  second_nearest_neighbor_list_.emplace_back(index);
}
void Atom::AppendThirdNearestNeighborList(int index) {
  third_nearest_neighbor_list_.emplace_back(index);
}
const std::vector<int> &Atom::GetSecondNearestNeighborList() const {
  return second_nearest_neighbor_list_;
}
const std::vector<int> &Atom::GetThirdNearestNeighborList() const {
  return third_nearest_neighbor_list_;
}

Vector3 GetRelativeDistanceVector(const Atom &first, const Atom &second) {
  Vector3 relative_distance_vector = second.GetRelativePosition() - first.GetRelativePosition();
  // periodic boundary conditions
  for (const auto kDim : All_Dimensions) {
    if (relative_distance_vector[kDim] >= 0.5)
      relative_distance_vector[kDim] -= 1;
    else if (relative_distance_vector[kDim] < -0.5)
      relative_distance_vector[kDim] += 1;
  }
  return relative_distance_vector;
}

} // namespace kn