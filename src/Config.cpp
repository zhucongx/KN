#include"Config.h"

#include <iostream>
#include <utility>

namespace kn {
Config::Config() = default;
Config::Config(const Matrix33 &basis, int atom_size) : basis_(basis) {
  atom_list_.reserve(atom_size);
}

bool Config::operator<(const Config &rhs) const {
  return energy_ < rhs.energy_;
}

void Config::ConvertRelativeToCartesian() {
  for (auto &atom : atom_list_) {
    atom.SetCartesianPosition(atom.GetRelativePosition() * basis_);
  }
}

void Config::ConvertCartesianToRelative() {
  auto inverse_basis = InverseMatrix33(basis_);
  for (auto &atom : atom_list_) {
    atom.SetRelativePosition(atom.GetCartesianPosition() * inverse_basis);
  }
}

void Config::UpdateNeighbors(double first_r_cutoff, double second_r_cutoff, double third_r_cutoff) {
  if (neighbor_found_)
    return;

  double first_r_cutoff_square = first_r_cutoff * first_r_cutoff;
  double second_r_cutoff_square = second_r_cutoff * second_r_cutoff;
  double third_r_cutoff_square = third_r_cutoff * third_r_cutoff;

  for (auto it1 = atom_list_.begin(); it1 < atom_list_.end(); ++it1) {
    for (auto it2 = atom_list_.begin(); it2 < it1; ++it2) {
      Vector3 absolute_distance_vector = GetRelativeDistanceVector(*it1, *it2) * basis_;
      if (absolute_distance_vector[kXDimension] > third_r_cutoff_square)
        continue;
      if (absolute_distance_vector[kYDimension] > third_r_cutoff_square)
        continue;
      if (absolute_distance_vector[kZDimension] > third_r_cutoff_square)
        continue;
      double absolute_distance_square = Inner(absolute_distance_vector);
      if (absolute_distance_square <= third_r_cutoff_square){
        if (absolute_distance_square <= second_r_cutoff_square) {
          if (absolute_distance_square <= first_r_cutoff_square) {
            it1->AppendFirstNearestNeighborList(it2->GetId());
            it2->AppendFirstNearestNeighborList(it1->GetId());
          } else {
            it1->AppendSecondNearestNeighborList(it2->GetId());
            it2->AppendSecondNearestNeighborList(it1->GetId());
          }
        } else{
          it1->AppendThirdNearestNeighborList(it2->GetId());
          it2->AppendThirdNearestNeighborList(it1->GetId());
        }
      }
    }
  }
  neighbor_found_ = true;
}
void Config::WrapAtomRelative() {
  for (auto &atom : atom_list_) {
    atom.SetRelativePosition(atom.GetRelativePosition() - ElementFloor(atom.GetRelativePosition()));
    atom.SetCartesianPosition(atom.GetRelativePosition() * basis_);
  }
}
void Config::WrapAtomCartesian() {
  auto inverse_basis = InverseMatrix33(basis_);
  for (auto &atom : atom_list_) {
    atom.SetRelativePosition(atom.GetCartesianPosition() * inverse_basis);
    atom.SetRelativePosition(atom.GetRelativePosition() - ElementFloor(atom.GetRelativePosition()));
    atom.SetCartesianPosition(atom.GetRelativePosition() * basis_);
  }
}

// for better performance, shouldn't call Wrap function
void Config::MoveRelativeDistance(const Vector3 &distance_vector) {
  for (auto &atom : atom_list_) {
    atom.SetRelativePosition(atom.GetRelativePosition() + distance_vector);
    atom.SetRelativePosition(atom.GetRelativePosition() - ElementFloor(atom.GetRelativePosition()));
    atom.SetCartesianPosition(atom.GetRelativePosition() * basis_);
  }
}

void Config::MoveOneAtomRelativeDistance(int index,
                                         const Vector3 &distance_vector) {
  atom_list_[index].SetRelativePosition(atom_list_[index].GetRelativePosition() + distance_vector);
  atom_list_[index].SetRelativePosition(atom_list_[index].GetRelativePosition()
                                            - ElementFloor(atom_list_[index].GetRelativePosition()));
  atom_list_[index].SetCartesianPosition(atom_list_[index].GetRelativePosition() * basis_);
}

// void Config::MoveAbsoluteDistance(const Vector3 &distance_vector) {
//   for (auto &atom:atom_list_) {
//     atom.cartesian_position_ = atom.cartesian_position_ + distance_vector;
//   }
//   WrapAbsolutePosition();
// }

const double kMean = 0;
const double kStandardDeviation = 0.15;
const double kPerturbCutOff = 0.4;

void Config::Perturb(std::mt19937_64 &generator) {
  std::normal_distribution<double> distribution(kMean, kStandardDeviation);
  auto add_displacement = [&generator, &distribution](double &coordinate) {
    double displacement = distribution(generator);
    while (std::abs(displacement) > kPerturbCutOff) {
      displacement = distribution(generator);
    }
    coordinate += displacement;
  };
  for (auto &atom : atom_list_) {
    auto cartesian_position = atom.GetCartesianPosition();
    for (const auto kDim : All_Dimensions) {
      add_displacement(cartesian_position[kDim]);
    }
    atom.SetCartesianPosition(cartesian_position);
  }
  WrapAtomCartesian();
}

int Config::GetNumAtoms() const {
  return atom_list_.size();
}

const Matrix33 &Config::GetBasis() const {
  return basis_;
}

void Config::AppendAtom(const Atom &atom) {
  atom_list_.push_back(atom);
  element_list_map_[atom.GetType()].emplace_back(atom.GetId());
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

std::map<Bond, int> CountAllBonds(Config &config) {
  std::map<Bond, int> bonds_count_map;
  std::string type1, type2;
  auto atom_list = config.GetAtomList();
  for (const auto &atom : atom_list) {
    type1 = atom.GetType();
    for (const auto &atom2_id : atom.GetFirstNearestNeighborList()) {
      bonds_count_map[Bond{type1, atom_list[atom2_id].GetType()}]++;
    }
  }
  for (auto &bond_count : bonds_count_map) {
    bond_count.second /= 2;
  }
  return bonds_count_map;
}
} // namespace kn
