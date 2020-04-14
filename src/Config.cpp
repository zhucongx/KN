#include"Config.h"

namespace box {

const double kMean = 0;
const double kStandardDeviation = 0.15;
const double kPerturbCutOff = 0.4;

Config::Config() = default;

bool Config::operator<(const Config &rhs) const {
  return energy_ < rhs.energy_;
}

void Config::Initialize() {
  first_bravais_vector_ = {0, 0, 0};
  second_bravais_vector_ = {0, 0, 0};
  third_bravais_vector_ = {0, 0, 0};
  num_atoms_ = 0;
  energy_ = 0.0;
  scale_ = 1.0;
  atom_list_.clear();
  element_list_set_.clear();
  neighborFound = false;
}
bool Config::IsCubic() const {
  return first_bravais_vector_.x == second_bravais_vector_.y &&
      second_bravais_vector_.y == third_bravais_vector_.z &&
      first_bravais_vector_.y == 0 &&
      first_bravais_vector_.z == 0 &&
      second_bravais_vector_.x == 0 &&
      second_bravais_vector_.z == 0 &&
      third_bravais_vector_.x == 0 &&
      third_bravais_vector_.y == 0;
}
void Config::ConvertRelativeToAbsolute() {
  for (auto &atom:atom_list_) {
    atom.absolute_position_ =
        LinearTransform(atom.relative_position_,
                                      first_bravais_vector_,
                                      second_bravais_vector_,
                                      third_bravais_vector_);

  }
}
void Config::ConvertAbsoluteToRelative() {
  arma::mat bravais_matrix =
      {{first_bravais_vector_.x, first_bravais_vector_.y,
        first_bravais_vector_.z},
       {second_bravais_vector_.x, second_bravais_vector_.y,
        second_bravais_vector_.z},
       {third_bravais_vector_.x, third_bravais_vector_.y,
        third_bravais_vector_.z}};
  arma::mat inverse_matrix = arma::inv(bravais_matrix);

  Vector3<double> first_inverse_vector = {inverse_matrix(kXDim, kXDim),
                                  inverse_matrix(kXDim, kYDim),
                                  inverse_matrix(kXDim, kZDim)};
  Vector3<double> second_inverse_vector = {inverse_matrix(kYDim, kXDim),
                                   inverse_matrix(kYDim, kYDim),
                                   inverse_matrix(kYDim, kZDim)};
  Vector3<double> third_inverse_vector = {inverse_matrix(kZDim, kXDim),
                                  inverse_matrix(kZDim, kYDim),
                                  inverse_matrix(kZDim, kZDim)};
  for (auto &atom:atom_list_) {
    atom.relative_position_ =
        LinearTransform(atom.absolute_position_,
                                      first_inverse_vector,
                                      second_inverse_vector,
                                      third_inverse_vector);
  }
}

void Config::Perturb() {
  auto seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::mt19937_64 generator(seed);
  std::normal_distribution<double> distribution(kMean, kStandardDeviation);
  auto add_displacement = [&generator, &distribution](double &coordinate) {
    double displacement = distribution(generator);
    while (abs(displacement) > kPerturbCutOff) {
      displacement = distribution(generator);
    }
    coordinate += displacement;
  };
  for (auto &atom:atom_list_) {
    add_displacement(atom.absolute_position_.x);
    add_displacement(atom.absolute_position_.y);
    add_displacement(atom.absolute_position_.z);
  }
  ConvertAbsoluteToRelative();
}
void Config::UpdateNeighbors(double first_r_cutoff, double second_r_cutoff) {
  double first_r_cutoff_square = first_r_cutoff * first_r_cutoff;
  double second_r_cutoff_square = second_r_cutoff * second_r_cutoff;
  for (auto it1 = atom_list_.begin(); it1 < atom_list_.end(); ++it1) {
    for (auto it2 = atom_list_.begin(); it2 < it1; ++it2) {
      Vector3<double> absolute_distance_vector =
          LinearTransform(GetRelativeDistanceVector(*it1, *it2),
                                        first_bravais_vector_,
                                        second_bravais_vector_,
                                        third_bravais_vector_);
      if (absolute_distance_vector.x > second_r_cutoff_square)
        continue;
      if (absolute_distance_vector.y > second_r_cutoff_square)
        continue;
      if (absolute_distance_vector.z > second_r_cutoff_square)
        continue;
      double absolute_distance_square =
          InnerProduct(absolute_distance_vector);
      if (absolute_distance_square <= second_r_cutoff_square) {
        if (absolute_distance_square <= first_r_cutoff_square) {
          it1->first_nearest_neighbor_list_.emplace_back(it2->GetId());
          it2->first_nearest_neighbor_list_.emplace_back(it1->GetId());
        }
        it1->second_nearest_neighbor_list_.emplace_back(it2->GetId());
        it2->second_nearest_neighbor_list_.emplace_back(it1->GetId());
      }
    }
  }
  neighborFound = true;
}

void Config::WrapPeriodicAtomPosition() {
  auto check_periodic = [](double &coordinate) {
    coordinate -= floor(coordinate);
  };
  for (auto &atom:atom_list_) {
    check_periodic(atom.relative_position_.x);
    check_periodic(atom.relative_position_.y);
    check_periodic(atom.relative_position_.z);
  }
  ConvertRelativeToAbsolute();
}

void Config::ShiftAtomToCentral(const Atom::Rank &id) {
  Vector3<double> criterionDistance =
      atom_list_[id].relative_position_ - Vector3<double>{0.5, 0.5, 0.5};
  for (auto &atom:atom_list_) {
    atom.relative_position_ = atom.relative_position_ - criterionDistance;
  }
  WrapPeriodicAtomPosition();
}

void Config::MoveAbsoluteDistance(const Vector3<double> &distance_vector) {
  for (auto &atom:atom_list_) {
    atom.absolute_position_ = atom.absolute_position_ + distance_vector;
  }
  ConvertAbsoluteToRelative();
  WrapPeriodicAtomPosition();
}

std::map<std::string, int> Config::CountAllBonds(double r_cutoff) {
  if (!neighborFound)
    UpdateNeighbors(r_cutoff, r_cutoff);

  std::map<std::string, int> bonds_count_map;
  std::string type1, type2, bond;
  for (const auto &atom:atom_list_) {
    type1 = atom.GetType();
    for (const auto &atom2_id:atom.first_nearest_neighbor_list_) {
      type2 = atom_list_[atom2_id].GetType();
      if (type1.compare(type2) < 0) {
        bond = type1;
        bond += "-";
        bond += type2;
      } else {
        bond = type2;
        bond += "-";
        bond += type1;
      }
      bonds_count_map[bond]++;
    }
  }
  for (auto &[bond, count]:bonds_count_map) {
    count /= 2;
  }
  return bonds_count_map;
}

}// namespace box
