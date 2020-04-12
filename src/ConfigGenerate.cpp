#include "Config.h"

void Config::GenerateFCC(const double &lattice_constant_a,
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
                                z_reference / z_length);
        atom_list_.emplace_back(num_atoms_++, mass, element,
                                (x_reference + 0.5) / x_length,
                                y_reference / y_length,
                                (z_reference + 0.5) / z_length);
        atom_list_.emplace_back(num_atoms_++, mass, element,
                                x_reference / x_length,
                                (y_reference + 0.5) / y_length,
                                (z_reference + 0.5) / z_length);
      }
    }
  }
  ConvertRelativeToAbsolute();
}
void Config::GenerateBCC(const double &lattice_constant_a,
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

void Config::GenerateHCP(const double &lattice_constant_a,
                         const double &lattice_constant_c,
                         const std::string &element,
                         const Int3 &factors) {
  Initialize();
  double mass = elem_info::FindMass(element);
  box_.SetFirstBravaisVector({lattice_constant_a * factors[kXDim], 0, 0});
  box_.SetSecondBravaisVector({-0.5 * lattice_constant_a * factors[kYDim],
                               0.5 * sqrt(3) * lattice_constant_a
                                   * factors[kYDim], 0});
  box_.SetThirdBravaisVector({0, 0, lattice_constant_c * factors[kZDim]});
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
                                (x_reference + 1.0 / 3.0) / x_length,
                                (y_reference + 2.0 / 3.0) / y_length,
                                (z_reference + 0.5) / z_length);
      }
    }
  }
  ConvertRelativeToAbsolute();
}
// void Config::embedCluster(const std::pair<std::string, std::string> &Elems,
//                           const FCCEmbededCluster::occupInfo_256 &o256,
//                           const int &i) {
//   for (int j = 0; j < o256.mapping[i].size(); ++j) {
//     if (o256.mapping[i][j] == 1) {
//       atoms[j].setType(Elems.first) ;
//     } else if (o256.mapping[i][j] == 2) {
//       atoms[j].setType(Elems.second);
//     }
//   }
// }