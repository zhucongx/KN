#include"Config.h"

#include <iostream>
#include <utility>

#include "AtomUtility.h"

namespace kn {

Config::Config() = default;

bool Config::operator<(const Config &rhs) const {
  return energy_ < rhs.energy_;
}

void Config::ConvertRelativeToCartesian() {
  for (auto &atom : atom_list_) {
    atom.cartesian_position_ = atom.relative_position_ * basis_;
  }
}

void Config::ConvertCartesianToRelative() {
  auto inverse_basis = InverseMatrix33(basis_);
  for (auto &atom : atom_list_) {
    atom.relative_position_ = atom.cartesian_position_ * inverse_basis;
  }
}



void Config::UpdateNeighbors(double first_r_cutoff, double second_r_cutoff) {
  if (neighbor_found_)
    return;

  double first_r_cutoff_square = first_r_cutoff * first_r_cutoff;
  double second_r_cutoff_square = second_r_cutoff * second_r_cutoff;
  for (auto it1 = atom_list_.begin(); it1 < atom_list_.end(); ++it1) {
    for (auto it2 = atom_list_.begin(); it2 < it1; ++it2) {
      Vector3 absolute_distance_vector =
          AtomUtility::GetRelativeDistanceVector(*it1, *it2) * basis_;
      if (absolute_distance_vector[kXDimension] > second_r_cutoff_square)
        continue;
      if (absolute_distance_vector[kYDimension] > second_r_cutoff_square)
        continue;
      if (absolute_distance_vector[kZDimension] > second_r_cutoff_square)
        continue;
      double absolute_distance_square = Inner(absolute_distance_vector);
      if (absolute_distance_square <= second_r_cutoff_square) {
        if (absolute_distance_square <= first_r_cutoff_square) {
          it1->first_nearest_neighbor_list_.emplace_back(it2->id_);
          it2->first_nearest_neighbor_list_.emplace_back(it1->id_);
        } else {
          it1->second_nearest_neighbor_list_.emplace_back(it2->id_);
          it2->second_nearest_neighbor_list_.emplace_back(it1->id_);
        }
      }
    }
  }
  neighbor_found_ = true;
}

// for better performance, shouldn't call Wrap function
void Config::MoveRelativeDistance(const Vector3 &distance_vector) {
  for (auto &atom : atom_list_) {
    atom.relative_position_ += distance_vector;

    atom.relative_position_ -= ElementFloor(atom.relative_position_);

    atom.cartesian_position_ = atom.relative_position_ * basis_;
  }
}

void Config::MoveOneAtomRelativeDistance(const int &index,
                                         const Vector3 &distance_vector) {
  atom_list_[index].relative_position_ += distance_vector;
  atom_list_[index].relative_position_ -= ElementFloor(atom_list_[index].relative_position_);

  atom_list_[index].cartesian_position_ = atom_list_[index].relative_position_ * basis_;
}

// void Config::MoveAbsoluteDistance(const Vector3 &distance_vector) {
//   for (auto &atom:atom_list_) {
//     atom.cartesian_position_ = atom.cartesian_position_ + distance_vector;
//   }
//   WrapAbsolutePosition();
// }


int Config::GetNumAtoms() const {
  return atom_list_.size();
}

void Config::SetScale(double scale) {
  scale_ = scale;
}

double Config::GetScale() const {
  return scale_;
}

const Matrix33 &Config::GetBasis() const {
  return basis_;
}

void Config::SetBasis(const Matrix33 &basis) {
  basis_ = basis;
}

void Config::AppendAtom(const Atom &atom) {
  atom_list_.push_back(atom);
  element_list_map_[atom.type_].emplace_back(atom.id_);
}

const std::vector<Atom> &Config::GetAtomList() const {
  return atom_list_;
}

const std::map<std::string, std::vector<int>> &Config::GetElementListMap() const {
  return element_list_map_;
}

bool Config::IsNeighborFound() const {
  return neighbor_found_;
}

void Config::SetNeighborFound(bool neighbor_found) {
  neighbor_found_ = neighbor_found;
}

} // namespace kn
