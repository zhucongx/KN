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
  box_.Initialize();
  atom_list_.clear();
  vacancy_list_.clear();
}

void Config::ConvertRelativeToAbsolute() {
  auto first_bravais_vector = box_.GetFirstBravaisVector();
  auto second_bravais_vector = box_.GetSecondBravaisVector();
  auto third_bravais_vector = box_.GetThirdBravaisVector();

  for (auto &atom:atom_list_) {

    auto relative_position = atom.GetRelativePosition();
    Double3 absolute_position =
        double3_calc::LinearTransform(relative_position,
                                      first_bravais_vector,
                                      second_bravais_vector,
                                      third_bravais_vector);
    atom.SetAbsolutePosition(absolute_position);
  }
}
void Config::ConvertAbsoluteToRelative() {
  auto first_bravais_vector = box_.GetFirstBravaisVector();
  auto second_bravais_vector = box_.GetSecondBravaisVector();
  auto third_bravais_vector = box_.GetThirdBravaisVector();

  arma::mat bravais_matrix =
      {{first_bravais_vector.x, first_bravais_vector.y,
        first_bravais_vector.z},
       {second_bravais_vector.x, second_bravais_vector.y,
        second_bravais_vector.z},
       {third_bravais_vector.x, third_bravais_vector.y,
        third_bravais_vector.z}};
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
    Double3 absolute_position =
        atom.GetAbsolutePosition();
    Double3 relative_position =
        double3_calc::LinearTransform(absolute_position,
                                      first_inverse_vector,
                                      second_inverse_vector,
                                      third_inverse_vector);
    atom.SetRelativePosition(relative_position);
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
    auto absolute_position = atom.GetAbsolutePosition();
    add_displacement(absolute_position.x);
    add_displacement(absolute_position.y);
    add_displacement(absolute_position.z);
    atom.SetAbsolutePosition(absolute_position);
  }
  ConvertAbsoluteToRelative();
}
void Config::UpdateNeighbors(double first_r_cutoff, double second_r_cutoff) {
  auto first_bravais_vector = box_.GetFirstBravaisVector();
  auto second_bravais_vector = box_.GetSecondBravaisVector();
  auto third_bravais_vector = box_.GetThirdBravaisVector();
  double first_r_cutoff_square = first_r_cutoff * first_r_cutoff;
  double second_r_cutoff_square = second_r_cutoff * second_r_cutoff;
  for (int i = 0; i < num_atoms_; i++) {
    for (int j = 0; j < i; j++) {
      Double3 absolute_distance_vector =
          double3_calc::LinearTransform(GetRelativeDistanceVector(i, j),
                                        first_bravais_vector,
                                        second_bravais_vector,
                                        third_bravais_vector);
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
          atom_list_[i].first_nearest_neighbor_list_.emplace_back(j);
          atom_list_[j].first_nearest_neighbor_list_.emplace_back(i);
        }
        atom_list_[i].second_nearest_neighbor_list_.emplace_back(j);
        atom_list_[j].second_nearest_neighbor_list_.emplace_back(i);
      }
    }
  }
}

Double3 Config::GetRelativeDistanceVector(int first, int second) const {
  return box::GetRelativeDistanceVector(atom_list_[first], atom_list_[second]);
}

}// namespace box
