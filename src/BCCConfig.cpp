#include "BCCConfig.h"
void BCCConfig::GenerateBCC(const double &lattice_constant_a,
                         const std::string &element,
                         const Int3 &factors) {
  Initialize();
  double mass = elem_info::FindMass(element);
  box_.SetFirstBravaisVector({lattice_constant_a * factors[kXDim], 0, 0});
  box_.SetSecondBravaisVector({0, lattice_constant_a * factors[kYDim], 0});
  box_.SetThirdBravaisVector({0, 0, lattice_constant_a * factors[kZDim]});
  box_.SetScale(1.0);
  num_atoms_ = 0;
  auto x_length = static_cast<double>(factors[kXDim]);
  auto y_length = static_cast<double>(factors[kYDim]);
  auto z_length = static_cast<double>(factors[kZDim]);
  for (int k = 0; k < factors[kZDim]; ++k) {
    for (int j = 0; j < factors[kYDim]; ++j) {
      for (int i = 0; i < factors[kXDim]; ++i) {
        auto x_reference = static_cast<double>(i);
        auto y_reference = static_cast<double>(j);
        auto z_reference = static_cast<double>(k);
        atom_list_.emplace_back(num_atoms_++, mass, element,
                                x_reference / x_length,
                                y_reference / y_length,
                                z_reference / z_length);
        atom_list_.emplace_back(num_atoms_++, mass, element,
                                (x_reference + 0.5) / x_length,
                                (y_reference + 0.5) / y_length,
                                (z_reference + 0.5) / z_length);
      }
    }
  }
  ConvertRelativeToAbsolute();
}
void BCCConfig::UpdateNeighbors(double first_r_cutoff, double second_r_cutoff) {
  Config::UpdateNeighbors(first_r_cutoff, second_r_cutoff);
}
