#include "Atom.h"
namespace cfg {
Atom::Atom() = default;

Atom::Atom(int id, double mass, std::string type, double x, double y, double z) :
    id_(id),
    mass_(mass),
    type_(std::move(type)),
    cartesian_position_{x, y, z},
    relative_position_{x, y, z} {
  first_nearest_neighbors_list_.reserve(Al_const::kNumFirstNearestNeighbors);
  second_nearest_neighbors_list_.reserve(Al_const::kNumSecondNearestNeighbors);
  third_nearest_neighbors_list_.reserve(Al_const::kNumThirdNearestNeighbors);
}
Atom::Atom(int id, double mass, std::string type, Vector3 position) :
    id_(id),
    mass_(mass),
    type_(std::move(type)),
    cartesian_position_(position),
    relative_position_(position) {
  first_nearest_neighbors_list_.reserve(Al_const::kNumFirstNearestNeighbors);
  second_nearest_neighbors_list_.reserve(Al_const::kNumSecondNearestNeighbors);
  third_nearest_neighbors_list_.reserve(Al_const::kNumThirdNearestNeighbors);
}

bool operator==(const Atom &lhs, const Atom &rhs) {
  return lhs.id_ == rhs.id_;
}
bool operator<(const Atom &lhs, const Atom &rhs) {
  return lhs.id_ < rhs.id_;
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
int Atom::GetId() const {
  return id_;
}
void Atom::SetId(int id) {
  id_ = id;
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
void Atom::AppendFirstNearestNeighborsList(int index) {
  first_nearest_neighbors_list_.emplace_back(index);
}
void Atom::AppendSecondNearestNeighborsList(int index) {
  second_nearest_neighbors_list_.emplace_back(index);
}
void Atom::AppendThirdNearestNeighborsList(int index) {
  third_nearest_neighbors_list_.emplace_back(index);
}

const std::vector<int> &Atom::GetFirstNearestNeighborsList() const {
  return first_nearest_neighbors_list_;
}
const std::vector<int> &Atom::GetSecondNearestNeighborsList() const {
  return second_nearest_neighbors_list_;
}
const std::vector<int> &Atom::GetThirdNearestNeighborsList() const {
  return third_nearest_neighbors_list_;
}
std::unordered_set<int> Atom::GetFirstAndSecondNeighborsSet() const {
  std::unordered_set<int> near_neighbors_hashset;
  std::copy(first_nearest_neighbors_list_.begin(),
            first_nearest_neighbors_list_.end(),
            std::inserter(near_neighbors_hashset, near_neighbors_hashset.end()));
  std::copy(second_nearest_neighbors_list_.begin(),
            second_nearest_neighbors_list_.end(),
            std::inserter(near_neighbors_hashset, near_neighbors_hashset.end()));
  std::copy(third_nearest_neighbors_list_.begin(),
            third_nearest_neighbors_list_.end(),
            std::inserter(near_neighbors_hashset, near_neighbors_hashset.end()));
  return near_neighbors_hashset;
}
std::unordered_set<int> Atom::GetFirstAndSecondThirdNeighborsSet() const {
  std::unordered_set<int> near_neighbors_hashset;
  std::copy(first_nearest_neighbors_list_.begin(),
            first_nearest_neighbors_list_.end(),
            std::inserter(near_neighbors_hashset, near_neighbors_hashset.end()));
  std::copy(second_nearest_neighbors_list_.begin(),
            second_nearest_neighbors_list_.end(),
            std::inserter(near_neighbors_hashset, near_neighbors_hashset.end()));
  std::copy(third_nearest_neighbors_list_.begin(),
            third_nearest_neighbors_list_.end(),
            std::inserter(near_neighbors_hashset, near_neighbors_hashset.end()));
  return near_neighbors_hashset;
}

void Atom::CleanNeighborsLists() {
  first_nearest_neighbors_list_.clear();
  second_nearest_neighbors_list_.clear();
  third_nearest_neighbors_list_.clear();
}
void AtomsJump(Atom &lhs, Atom &rhs) {
  std::swap(lhs.relative_position_, rhs.relative_position_);
  std::swap(lhs.cartesian_position_, rhs.cartesian_position_);

  std::swap(lhs.first_nearest_neighbors_list_, rhs.first_nearest_neighbors_list_);
  std::replace(lhs.first_nearest_neighbors_list_.begin(),
               lhs.first_nearest_neighbors_list_.end(),
               lhs.GetId(), rhs.GetId());
  std::replace(rhs.first_nearest_neighbors_list_.begin(),
               rhs.first_nearest_neighbors_list_.end(),
               rhs.GetId(), lhs.GetId());

  std::swap(lhs.second_nearest_neighbors_list_, rhs.second_nearest_neighbors_list_);
  std::replace(lhs.second_nearest_neighbors_list_.begin(),
               lhs.second_nearest_neighbors_list_.end(),
               lhs.GetId(), rhs.GetId());
  std::replace(rhs.second_nearest_neighbors_list_.begin(),
               rhs.second_nearest_neighbors_list_.end(),
               rhs.GetId(), lhs.GetId());

  std::swap(lhs.third_nearest_neighbors_list_, rhs.third_nearest_neighbors_list_);
  std::replace(lhs.third_nearest_neighbors_list_.begin(),
               lhs.third_nearest_neighbors_list_.end(),
               lhs.GetId(), rhs.GetId());
  std::replace(rhs.third_nearest_neighbors_list_.begin(),
               rhs.third_nearest_neighbors_list_.end(),
               rhs.GetId(), lhs.GetId());
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

} // namespace cfg