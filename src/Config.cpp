#include"Config.h"

const double kMean = 0;
const double kStandardDeviation = 0.15;
const double kCutOff = 0.4;

Config::Config() = default;


void Config::Initialize() {
  num_atoms_ = 0;
  energy_ = 0.0;
  box_.Initialize();
  atom_list_.clear();
  vacancy_list_.clear();
}

bool Config::operator<(const Config &rhs) const {
  return energy_ < rhs.energy_;
}


void Config::ConvertRelativeToAbsolute() {
  auto first_bravais_vector = box_.GetFirstBravaisVector();
  auto second_bravais_vector = box_.GetSecondBravaisVector();
  auto third_bravais_vector = box_.GetThirdBravaisVector();

  for (auto &atom:atom_list_) {
    std::array<double, kDimension> absolute_position{};
    auto relative_position = atom.GetRelativePosition();
    for (const auto &i : {kXDim, kYDim, kZDim}) {
      absolute_position[i] =
          relative_position[kXDim] * first_bravais_vector[i]
              + relative_position[kYDim] * second_bravais_vector[i]
              + relative_position[kZDim] * third_bravais_vector[i];
    }
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
  for (auto &atom:atom_list_) {
    std::array<double, kDimension> absolute_position =
        atom.GetAbsolutePosition();
    std::array<double, kDimension> relative_position{};
    arma::vec b = {absolute_position[kXDim],
                   absolute_position[kYDim],
                   absolute_position[kZDim]};
    arma::vec x = solve(bravais_matrix, b);
    for (const auto &i : {kXDim, kYDim, kZDim}) {
      relative_position[i] = x[i];
    }
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