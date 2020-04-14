#include "Config.h"

namespace box {

void Config::GenerateFCC(const double &lattice_constant_a,
                         const std::string &element,
                         const Vector3<int> &factors) {
  Initialize();
  double mass = elem_info::FindMass(element);
  first_bravais_vector_ = {lattice_constant_a * factors.x, 0, 0};
  second_bravais_vector_ = {0, lattice_constant_a * factors.y, 0};
  third_bravais_vector_ = {0, 0, lattice_constant_a * factors.z};
  scale_ = 1.0;
  num_atoms_ = 0;
  auto x_length = static_cast<double>(factors.x);
  auto y_length = static_cast<double>(factors.y);
  auto z_length = static_cast<double>(factors.z);
  for (int k = 0; k < factors.z; ++k) {
    for (int j = 0; j < factors.y; ++j) {
      for (int i = 0; i < factors.x; ++i) {
        auto x_reference = static_cast<double>(i);
        auto y_reference = static_cast<double>(j);
        auto z_reference = static_cast<double>(k);
        atom_list_.emplace_back(num_atoms_, mass, element,
                                x_reference / x_length,
                                y_reference / y_length,
                                z_reference / z_length);
        element_list_set_[element].emplace_back(num_atoms_++);
        atom_list_.emplace_back(num_atoms_, mass, element,
                                (x_reference + 0.5) / x_length,
                                (y_reference + 0.5) / y_length,
                                z_reference / z_length);
        element_list_set_[element].emplace_back(num_atoms_++);
        atom_list_.emplace_back(num_atoms_, mass, element,
                                (x_reference + 0.5) / x_length,
                                y_reference / y_length,
                                (z_reference + 0.5) / z_length);
        element_list_set_[element].emplace_back(num_atoms_++);
        atom_list_.emplace_back(num_atoms_, mass, element,
                                x_reference / x_length,
                                (y_reference + 0.5) / y_length,
                                (z_reference + 0.5) / z_length);
        element_list_set_[element].emplace_back(num_atoms_++);
      }
    }
  }
  ConvertRelativeToAbsolute();
}
void Config::GenerateHCP(const double &lattice_constant_a,
                         const double &lattice_constant_c,
                         const std::string &element,
                         const Vector3<int> &factors) {
  Initialize();
  double mass = elem_info::FindMass(element);
  first_bravais_vector_ = {lattice_constant_a * factors.x, 0, 0};
  second_bravais_vector_ = {-0.5 * lattice_constant_a * factors.y,
                            0.5 * sqrt(3) * lattice_constant_a * factors.y,
                            0};
  third_bravais_vector_ = {0, 0, lattice_constant_c * factors.z};
  scale_ = 1.0;
  num_atoms_ = 0;

  auto x_length = static_cast<double>(factors.x);
  auto y_length = static_cast<double>(factors.y);
  auto z_length = static_cast<double>(factors.z);
  for (int k = 0; k < factors.z; ++k) {
    for (int j = 0; j < factors.y; ++j) {
      for (int i = 0; i < factors.x; ++i) {
        auto x_reference = static_cast<double>(i);
        auto y_reference = static_cast<double>(j);
        auto z_reference = static_cast<double>(k);
        atom_list_.emplace_back(num_atoms_, mass, element,
                                x_reference / x_length,
                                y_reference / y_length,
                                z_reference / z_length);
        element_list_set_[element].emplace_back(num_atoms_++);
        atom_list_.emplace_back(num_atoms_, mass, element,
                                (x_reference + 1.0 / 3.0) / x_length,
                                (y_reference + 2.0 / 3.0) / y_length,
                                (z_reference + 0.5) / z_length);
        element_list_set_[element].emplace_back(num_atoms_++);
      }
    }
  }
  ConvertRelativeToAbsolute();
}

}// namespace box
