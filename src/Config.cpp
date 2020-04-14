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
  num_atoms_ = 0;
  energy_ = 0.0;
  scale_ = 1.0;
  first_bravais_vector_ = {0, 0, 0};
  second_bravais_vector_ = {0, 0, 0};
  third_bravais_vector_ = {0, 0, 0};
  atom_list_.clear();
  vacancy_list_.clear();
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
        double3_calc::LinearTransform(atom.relative_position_,
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

  Double3 first_inverse_vector = {inverse_matrix(kXDim, kXDim),
                                  inverse_matrix(kXDim, kYDim),
                                  inverse_matrix(kXDim, kZDim)};
  Double3 second_inverse_vector = {inverse_matrix(kYDim, kXDim),
                                   inverse_matrix(kYDim, kYDim),
                                   inverse_matrix(kYDim, kZDim)};
  Double3 third_inverse_vector = {inverse_matrix(kZDim, kXDim),
                                  inverse_matrix(kZDim, kYDim),
                                  inverse_matrix(kZDim, kZDim)};
  for (auto &atom:atom_list_) {
    atom.relative_position_ =
        double3_calc::LinearTransform(atom.absolute_position_,
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
  for (auto &atom : atom_list_) {
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
      Double3 absolute_distance_vector =
          double3_calc::LinearTransform(GetRelativeDistanceVector(*it1, *it2),
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
          double3_calc::InnerProduct(absolute_distance_vector);
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
}

void Config::ShiftAtomToCentral(const Rank &id) {
  Double3 criterionDistance = atom_list_[id].relative_position_;
  for (auto &atom:atom_list_) {
    atom.relative_position_ = double3_calc::Subtract(atom.relative_position_,
                                                     criterionDistance);
  }
}

}// namespace box
