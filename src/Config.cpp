#include"Config.h"

#include <iostream>
#include <utility>
#include <random>
#include <chrono>

namespace kn {

const double kMean = 0;
const double kStandardDeviation = 0.15;
const double kPerturbCutOff = 0.4;

Config::Config() = default;

bool Config::operator<(const Config &rhs) const {
  return energy_ < rhs.energy_;
}

void Config::Initialize() {
  basis_ = {{{0, 0, 0},
             {0, 0, 0},
             {0, 0, 0}}};
  energy_ = 0.0;
  scale_ = 1.0;
  atom_list_.clear();
  element_list_map_.clear();
  neighbor_found_ = false;
}

bool Config::IsCubic() const {
  return basis_[kXDimension][kXDimension] == basis_[kYDimension][kYDimension] &&
      basis_[kYDimension][kYDimension] == basis_[kZDimension][kZDimension] &&
      basis_[kXDimension][kYDimension] == 0 &&
      basis_[kXDimension][kZDimension] == 0 &&
      basis_[kYDimension][kXDimension] == 0 &&
      basis_[kYDimension][kZDimension] == 0 &&
      basis_[kZDimension][kXDimension] == 0 &&
      basis_[kZDimension][kYDimension] == 0;
}

void Config::ConvertRelativeToCartesian() {
  for (auto &atom:atom_list_) {
    atom.cartesian_position_ = atom.relative_position_ * basis_;
  }
}

void Config::ConvertCartesianToRelative() {
  auto inverse_basis = InverseMatrix33(basis_);
  for (auto &atom:atom_list_) {
    atom.relative_position_ = atom.cartesian_position_ * inverse_basis;
  }
}

void Config::Perturb() {
  auto seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::mt19937_64 generator(seed);
  std::normal_distribution<double> distribution(kMean, kStandardDeviation);
  auto add_displacement = [&generator, &distribution](double &coordinate) {
    double displacement = distribution(generator);
    while (std::abs(displacement) > kPerturbCutOff) {
      displacement = distribution(generator);
    }
    coordinate += displacement;
  };
  for (auto &atom:atom_list_) {
    add_displacement(atom.cartesian_position_[kXDimension]);
    add_displacement(atom.cartesian_position_[kYDimension]);
    add_displacement(atom.cartesian_position_[kZDimension]);
  }
  ConvertCartesianToRelative();
}

void Config::UpdateNeighbors(double first_r_cutoff, double second_r_cutoff) {
  if (neighbor_found_)
    return;

  double first_r_cutoff_square = first_r_cutoff * first_r_cutoff;
  double second_r_cutoff_square = second_r_cutoff * second_r_cutoff;
  for (auto it1 = atom_list_.begin(); it1 < atom_list_.end(); ++it1) {
    for (auto it2 = atom_list_.begin(); it2 < it1; ++it2) {
      Vector3 absolute_distance_vector =
          GetRelativeDistanceVector(*it1, *it2) * basis_;
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
        }
        it1->second_nearest_neighbor_list_.emplace_back(it2->id_);
        it2->second_nearest_neighbor_list_.emplace_back(it1->id_);
      }
    }
  }
  neighbor_found_ = true;
}

// void Config::WrapRelativePosition() {
//   for (auto &atom:atom_list_) {
//     atom.relative_position_ -= ElementFloor(atom.relative_position_);
//
//     atom.cartesian_position_ = atom.relative_position_ * basis_;
//   }
// }

// void Config::WrapAbsolutePosition() {
//   ConvertCartesianToRelative();
//   WrapRelativePosition();
// }

// void Config::ShiftAtomToCentral(const Atom::Rank &id) {
//   Vector3 criterionDistance =
//       atom_list_[id].relative_position_ - Vector3{0.5, 0.5, 0.5};
//   for (auto &atom:atom_list_) {
//     atom.relative_position_ = atom.relative_position_ - criterionDistance;
//   }
//   WrapRelativePosition();
// }
// for better performance, shouldn't call Wrap function
void Config::MoveRelativeDistance(const Vector3 &distance_vector) {
  for (auto &atom:atom_list_) {
    atom.relative_position_ += distance_vector;

    atom.relative_position_ -= ElementFloor(atom.relative_position_);

    atom.cartesian_position_ = atom.relative_position_ * basis_;
  }
}
void Config::MoveOneAtomRelativeDistance(const Atom::Rank &index,
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

std::map<Bond, int> Config::CountAllBonds(double r_cutoff) {
  if (!neighbor_found_)
    UpdateNeighbors(r_cutoff, r_cutoff);

  std::map<Bond, int> bonds_count_map;
  std::string type1, type2;
  for (const auto &atom:atom_list_) {
    type1 = atom.type_;
    for (const auto &atom2_id:atom.first_nearest_neighbor_list_) {
      bonds_count_map[Bond{type1, atom_list_[atom2_id].type_}]++;
    }
  }
  for (auto &bond_count:bonds_count_map) {
    bond_count.second /= 2;
  }
  return bonds_count_map;
}

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
const std::map<std::string, std::vector<Atom::Rank>> &Config::GetElementListMap() const {
  return element_list_map_;
}


}// namespace kn
