#include "ConfigUtility.h"
namespace kn {

const double kMean = 0;
const double kStandardDeviation = 0.15;
const double kPerturbCutOff = 0.4;

std::map<Bond, int> ConfigUtility::CountAllBonds(Config &config, double r_cutoff) {
  if (!config.IsNeighborFound())
    config.UpdateNeighbors(r_cutoff, r_cutoff);

  std::map<Bond, int> bonds_count_map;
  std::string type1, type2;
  auto atom_list = config.GetAtomList();
  for (const auto &atom : atom_list) {
    type1 = atom.type_;
    for (const auto &atom2_id : atom.first_nearest_neighbor_list_) {
      bonds_count_map[Bond{type1, atom_list[atom2_id].type_}]++;
    }
  }
  for (auto &bond_count : bonds_count_map) {
    bond_count.second /= 2;
  }
  return bonds_count_map;
}
void ConfigUtility::Perturb(Config &config, std::mt19937_64 &generator) {
  std::normal_distribution<double> distribution(kMean, kStandardDeviation);
  auto add_displacement = [&generator, &distribution](double &coordinate) {
    double displacement = distribution(generator);
    while (std::abs(displacement) > kPerturbCutOff) {
      displacement = distribution(generator);
    }
    coordinate += displacement;
  };
  for (auto &atom : config.atom_list_) {
    add_displacement(atom.cartesian_position_[kXDimension]);
    add_displacement(atom.cartesian_position_[kYDimension]);
    add_displacement(atom.cartesian_position_[kZDimension]);
  }
  config.ConvertCartesianToRelative();
}

} // namespace kn
