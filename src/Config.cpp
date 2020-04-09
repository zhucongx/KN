#include"Config.h"

constexpr double kMean = 0;
constexpr double kStandardDeviation = 0.15;
constexpr double kCutOff = 0.4;

Config::Config() = default;
Config::~Config() = default;

void Config::Initialize() {
  num_atoms_ = 0;
  scale_ = 1.0;
  energy_ = 0.0;
  first_bravais_vector_.fill(0.0);
  second_bravais_vector_.fill(0.0);
  third_bravais_vector_.fill(0.0);
  atom_list_.clear();
  vacancy_list_.clear();
}

bool Config::operator<(const Config &rhs) const {
  return energy_ < rhs.energy_;
}


void Config::ConvertRelativeToAbsolute() {
  for (auto &atom:atom_list_) {
    std::array<double, kDimension> absolute_position{};
    std::array<double, kDimension> relative_position =
        atom.GetRelativePosition();
    for (const auto &i : {kXDim, kYDim, kZDim}) {
      absolute_position[i] =
          relative_position[kXDim] * first_bravais_vector_[i]
              + relative_position[kYDim] * second_bravais_vector_[i]
              + relative_position[kZDim] * third_bravais_vector_[i];
    }
    atom.SetAbsolutePosition(absolute_position);
  }
}
void Config::ConvertAbsoluteToRelative() {
  arma::mat bm = {{first_bravais_vector_[kXDim], first_bravais_vector_[kYDim],
                   first_bravais_vector_[kZDim]},
                  {second_bravais_vector_[kXDim], second_bravais_vector_[kYDim],
                   second_bravais_vector_[kZDim]},
                  {third_bravais_vector_[kXDim], third_bravais_vector_[kYDim],
                   third_bravais_vector_[kZDim]}};
  for (auto &atom:atom_list_) {
    std::array<double, kDimension> absolute_position =
        atom.GetAbsolutePosition();
    std::array<double, kDimension> relative_position{};
    arma::vec b = {absolute_position[kXDim],
                   absolute_position[kYDim],
                   absolute_position[kZDim]};
    arma::vec x = solve(bm, b);
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