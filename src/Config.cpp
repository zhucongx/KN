#include"Config.h"
const double kMean = 0;
const double kStandardDeviation = 0.15;
const double kCutOff = 0.4;

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
      {{first_bravais_vector[kXDim], first_bravais_vector[kYDim],
        first_bravais_vector[kZDim]},
       {second_bravais_vector[kXDim], second_bravais_vector[kYDim],
        second_bravais_vector[kZDim]},
       {third_bravais_vector[kXDim], third_bravais_vector[kYDim],
        third_bravais_vector[kZDim]}};
  arma::mat inverse_bravais_matrix = arma::inv(bravais_matrix);
  Double3 first_inverse_bravais_vector, second_inverse_bravais_vector,
      third_inverse_bravais_vector;
  for (const auto &i : {kXDim, kYDim, kZDim}) {
    first_inverse_bravais_vector[i] = inverse_bravais_matrix(kXDim, i);
    second_inverse_bravais_vector[i] = inverse_bravais_matrix(kYDim, i);
    third_inverse_bravais_vector[i] = inverse_bravais_matrix(kZDim, i);
  }

  for (auto &atom:atom_list_) {
    Double3 absolute_position =
        atom.GetAbsolutePosition();
    Double3 relative_position =
        double3_calc::LinearTransform(absolute_position,
                                      first_inverse_bravais_vector,
                                      second_inverse_bravais_vector,
                                      third_inverse_bravais_vector);
    atom.SetRelativePosition(relative_position);
  }
}

void Config::Perturb() {
  auto seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::mt19937_64 generator(seed);
  std::normal_distribution<double> distribution(kMean, kStandardDeviation);
  for (auto &atom : atom_list_) {
    auto absolute_position = atom.GetAbsolutePosition();
    for (const auto &i : {kXDim, kYDim, kZDim}) {
      double displacement = distribution(generator);
      while (abs(displacement) > kCutOff) {
        displacement = distribution(generator);
      }
      absolute_position[i] += displacement;
    }
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
      if (absolute_distance_vector[kXDim] > second_r_cutoff_square)
        continue;
      if (absolute_distance_vector[kYDim] > second_r_cutoff_square)
        continue;
      if (absolute_distance_vector[kZDim] > second_r_cutoff_square)
        continue;
      double absolute_distance_square =
          double3_calc::InnerProduct(absolute_distance_vector);
      if (absolute_distance_square <= second_r_cutoff_square) {
        if (absolute_distance_square <= first_r_cutoff_square) {
          atom_list_[i].first_nearest_neighbor_list_.emplace_back(j);
          atom_list_[j].first_nearest_neighbor_list_.emplace_back(i);
        }
        atom_list_[i].second_near_neighbor_list_.emplace_back(j);
        atom_list_[j].second_near_neighbor_list_.emplace_back(i);
      }
    }
  }
}

Double3 Config::GetRelativeDistanceVector(int first, int second) const {
  auto atom1_relative_position = atom_list_[first].GetRelativePosition();
  auto atom2_relative_position = atom_list_[second].GetRelativePosition();
  Double3 relative_distance_vector =
      {atom2_relative_position[kXDim] - atom1_relative_position[kXDim],
       atom2_relative_position[kYDim] - atom1_relative_position[kYDim],
       atom2_relative_position[kZDim] - atom1_relative_position[kZDim]};

  // periodic boundary conditions
  for (auto &relative_distance_vector_element:relative_distance_vector) {
    if (relative_distance_vector_element >= 0.5)
      relative_distance_vector_element -= 1;
    else if (relative_distance_vector_element < -0.5)
      relative_distance_vector_element += 1;
  }
  return relative_distance_vector;
}


