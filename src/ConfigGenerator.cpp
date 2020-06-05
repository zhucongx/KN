#include "ConfigGenerator.h"

namespace kn {

Config ConfigGenerator::GenerateFCC(double lattice_constant_a,
                                    const std::string &element,
                                    const std::array<int, kDimension> &factors) {
  Config config;
  double mass = elem_info::FindMass(element);
  config.SetBasis({
      {
          {lattice_constant_a * factors[kXDimension], 0, 0},
          {0, lattice_constant_a * factors[kYDimension], 0},
          {0, 0, lattice_constant_a * factors[kZDimension]}
      }
  });
  config.SetScale(1.0);
  int atoms_counter = 0;
  auto x_length = static_cast<double>(factors[kXDimension]);
  auto y_length = static_cast<double>(factors[kYDimension]);
  auto z_length = static_cast<double>(factors[kZDimension]);
  for (int k = 0; k < factors[kZDimension]; ++k) {
    for (int j = 0; j < factors[kYDimension]; ++j) {
      for (int i = 0; i < factors[kXDimension]; ++i) {
        auto x_reference = static_cast<double>(i);
        auto y_reference = static_cast<double>(j);
        auto z_reference = static_cast<double>(k);
        config.AppendAtom({
            atoms_counter++, mass, element,
            x_reference / x_length,
            y_reference / y_length,
            z_reference / z_length
        });
        config.AppendAtom({
            atoms_counter++, mass, element,
            (x_reference + 0.5) / x_length,
            (y_reference + 0.5) / y_length,
            z_reference / z_length
        });
        config.AppendAtom({
            atoms_counter++, mass, element,
            (x_reference + 0.5) / x_length,
            y_reference / y_length,
            (z_reference + 0.5) / z_length
        });
        config.AppendAtom({
            atoms_counter++, mass, element,
            x_reference / x_length,
            (y_reference + 0.5) / y_length,
            (z_reference + 0.5) / z_length
        });
      }
    }
  }
  config.ConvertRelativeToCartesian();
  return config;
}

Config ConfigGenerator::GenerateBCC(double lattice_constant_a,
                                    const std::string &element,
                                    const std::array<int, kDimension> &factors) {
  Config config;
  double mass = elem_info::FindMass(element);
  config.SetBasis({
      {
          {lattice_constant_a * factors[kXDimension], 0, 0},
          {0, lattice_constant_a * factors[kYDimension], 0},
          {0, 0, lattice_constant_a * factors[kZDimension]}
      }
  });
  config.SetScale(1.0);
  int atoms_counter = 0;
  auto x_length = static_cast<double>(factors[kXDimension]);
  auto y_length = static_cast<double>(factors[kYDimension]);
  auto z_length = static_cast<double>(factors[kZDimension]);
  for (int k = 0; k < factors[kZDimension]; ++k) {
    for (int j = 0; j < factors[kYDimension]; ++j) {
      for (int i = 0; i < factors[kXDimension]; ++i) {
        auto x_reference = static_cast<double>(i);
        auto y_reference = static_cast<double>(j);
        auto z_reference = static_cast<double>(k);
        config.AppendAtom({
            atoms_counter++, mass, element,
            x_reference / x_length,
            y_reference / y_length,
            z_reference / z_length
        });
        config.AppendAtom({
            atoms_counter++, mass, element,
            (x_reference + 0.5) / x_length,
            (y_reference + 0.5) / y_length,
            (z_reference + 0.5) / z_length
        });
      }
    }
  }
  config.ConvertRelativeToCartesian();
  return config;
}

Config ConfigGenerator::GenerateHCP(double lattice_constant_a,
                                    double lattice_constant_c,
                                    const std::string &element,
                                    const std::array<int, kDimension> &factors) {
  Config config;
  double mass = elem_info::FindMass(element);
  config.SetBasis({
      {
          {lattice_constant_a * factors[kXDimension], 0, 0},
          {
              -0.5 * lattice_constant_a * factors[kYDimension],
              0.5 * sqrt(3) * lattice_constant_a * factors[kYDimension], 0
          },
          {0, 0, lattice_constant_c * factors[kZDimension]}
      }
  });
  config.SetScale(1.0);
  int atoms_counter = 0;

  auto x_length = static_cast<double>(factors[kXDimension]);
  auto y_length = static_cast<double>(factors[kYDimension]);
  auto z_length = static_cast<double>(factors[kZDimension]);
  for (int k = 0; k < factors[kZDimension]; ++k) {
    for (int j = 0; j < factors[kYDimension]; ++j) {
      for (int i = 0; i < factors[kXDimension]; ++i) {
        auto x_reference = static_cast<double>(i);
        auto y_reference = static_cast<double>(j);
        auto z_reference = static_cast<double>(k);
        config.AppendAtom({
            atoms_counter++, mass, element,
            x_reference / x_length,
            y_reference / y_length,
            z_reference / z_length
        });
        config.AppendAtom({
            atoms_counter, mass, element,
            (x_reference + 1.0 / 3.0) / x_length,
            (y_reference + 2.0 / 3.0) / y_length,
            (z_reference + 0.5) / z_length
        });
      }
    }
  }
  config.ConvertRelativeToCartesian();
  return config;
}

Config ConfigGenerator::GenerateL10(double lattice_constant_a,
                                    const std::vector<std::string> &element_list,
                                    const std::array<int, kDimension> &factors) {
  Matrix33 basis = {
      {
          {lattice_constant_a, 0, 0},
          {0, lattice_constant_a, 0},
          {0, 0, lattice_constant_a}
      }
  };
  std::vector<std::pair<std::string, Vector3>> element_position_list;
  element_position_list.emplace_back(element_list[0], Vector3{0.0, 0.0, 0.0});
  element_position_list.emplace_back(element_list[0], Vector3{0.5, 0.5, 0.0});
  element_position_list.emplace_back(element_list[1], Vector3{0.0, 0.5, 0.5});
  element_position_list.emplace_back(element_list[1], Vector3{0.5, 0.0, 0.5});

  return Duplicate(GenerateUnitCell(basis, element_position_list), factors);
}

Config ConfigGenerator::GenerateL12(double lattice_constant_a,
                                    const std::vector<std::string> &element_list,
                                    const std::array<int, kDimension> &factors) {
  Matrix33 basis = {
      {
          {lattice_constant_a, 0, 0},
          {0, lattice_constant_a, 0},
          {0, 0, lattice_constant_a}
      }
  };
  std::vector<std::pair<std::string, Vector3>> element_position_list;
  element_position_list.emplace_back(element_list[0], Vector3{0.0, 0.0, 0.0});
  element_position_list.emplace_back(element_list[1], Vector3{0.5, 0.5, 0.0});
  element_position_list.emplace_back(element_list[1], Vector3{0.0, 0.5, 0.5});
  element_position_list.emplace_back(element_list[1], Vector3{0.5, 0.0, 0.5});

  return Duplicate(GenerateUnitCell(basis, element_position_list), factors);
}

Config ConfigGenerator::GenerateL10star(double lattice_constant_a,
                                        const std::vector<std::string> &element_list,
                                        const std::array<int, kDimension> &factors) {
  Matrix33 basis = {
      {
          {lattice_constant_a * 2, 0, 0},
          {0, lattice_constant_a, 0},
          {0, 0, lattice_constant_a}
      }
  };
  std::vector<std::pair<std::string, Vector3>> element_position_list;
  element_position_list.emplace_back(element_list[0], Vector3{0.0, 0.0, 0.0});
  element_position_list.emplace_back(element_list[0], Vector3{0.5, 0.0, 0.0});
  element_position_list.emplace_back(element_list[1], Vector3{0.25, 0.5, 0.0});
  element_position_list.emplace_back(element_list[0], Vector3{0.75, 0.5, 0.0});
  element_position_list.emplace_back(element_list[0], Vector3{0.25, 0.0, 0.5});
  element_position_list.emplace_back(element_list[1], Vector3{0.75, 0.0, 0.5});
  element_position_list.emplace_back(element_list[1], Vector3{0.0, 0.5, 0.5});
  element_position_list.emplace_back(element_list[1], Vector3{0.5, 0.5, 0.5});

  return Duplicate(GenerateUnitCell(basis, element_position_list), factors);
}

Config ConfigGenerator::GenerateL12star(double lattice_constant_a,
                                        const std::vector<std::string> &element_list,
                                        const std::array<int, kDimension> &factors) {
  Matrix33 basis = {
      {
          {lattice_constant_a * 4, 0, 0},
          {0, lattice_constant_a, 0},
          {0, 0, lattice_constant_a}
      }
  };
  std::vector<std::pair<std::string, Vector3>> element_position_list;
  element_position_list.emplace_back(element_list[0], Vector3{0.0, 0.0, 0.0});
  element_position_list.emplace_back(element_list[1], Vector3{0.125, 0.5, 0.0});
  element_position_list.emplace_back(element_list[1], Vector3{0.25, 0.0, 0.0});
  element_position_list.emplace_back(element_list[1], Vector3{0.375, 0.5, 0.0});
  element_position_list.emplace_back(element_list[1], Vector3{0.5, 0.0, 0.0});
  element_position_list.emplace_back(element_list[1], Vector3{0.625, 0.5, 0.0});
  element_position_list.emplace_back(element_list[0], Vector3{0.75, 0.0, 0.0});
  element_position_list.emplace_back(element_list[1], Vector3{0.875, 0.5, 0.0});
  element_position_list.emplace_back(element_list[1], Vector3{0.0, 0.5, 0.5});
  element_position_list.emplace_back(element_list[1], Vector3{0.125, 0.0, 0.5});
  element_position_list.emplace_back(element_list[0], Vector3{0.25, 0.5, 0.5});
  element_position_list.emplace_back(element_list[1], Vector3{0.375, 0.0, 0.5});
  element_position_list.emplace_back(element_list[0], Vector3{0.5, 0.5, 0.5});
  element_position_list.emplace_back(element_list[1], Vector3{0.625, 0.0, 0.5});
  element_position_list.emplace_back(element_list[1], Vector3{0.75, 0.5, 0.5});
  element_position_list.emplace_back(element_list[1], Vector3{0.875, 0.0, 0.5});

  return Duplicate(GenerateUnitCell(basis, element_position_list), factors);
}

Config ConfigGenerator::GenerateZ1(double lattice_constant_a,
                                   const std::vector<std::string> &element_list,
                                   const std::array<int, kDimension> &factors) {
  Matrix33 basis = {
      {
          {lattice_constant_a, 0, 0},
          {0, lattice_constant_a, 0},
          {0, 0, lattice_constant_a * 2}
      }
  };
  std::vector<std::pair<std::string, Vector3>> element_position_list;
  element_position_list.emplace_back(element_list[0], Vector3{0.0, 0.0, 0.0});
  element_position_list.emplace_back(element_list[0], Vector3{0.5, 0.5, 0.0});
  element_position_list.emplace_back(element_list[1], Vector3{0.0, 0.5, 0.25});
  element_position_list.emplace_back(element_list[1], Vector3{0.5, 0.0, 0.25});
  element_position_list.emplace_back(element_list[1], Vector3{0.0, 0.0, 0.5});
  element_position_list.emplace_back(element_list[1], Vector3{0.5, 0.5, 0.5});
  element_position_list.emplace_back(element_list[1], Vector3{0.0, 0.5, 0.75});
  element_position_list.emplace_back(element_list[1], Vector3{0.5, 0.0, 0.75});

  return Duplicate(GenerateUnitCell(basis, element_position_list), factors);
}

Config ConfigGenerator::GenerateUnitCell(
    const Matrix33 &basis_matrix,
    const std::vector<std::pair<std::string, Vector3>> &type_position_list) {
  Config config;
  config.SetBasis(basis_matrix);
  config.SetScale(1.0);
  int atoms_counter = 0;
  for (const auto &[type, relative_position] : type_position_list) {
    config.AppendAtom({
        atoms_counter++,
        elem_info::FindMass(type),
        type,
        relative_position
    });
  }
  config.ConvertRelativeToCartesian();
  return config;
}

Config ConfigGenerator::Duplicate(const Config &in_config,
                                  const std::array<int, kDimension> &factors) {
  Config out_config;
  auto x_length = static_cast<double>(factors[kXDimension]);
  auto y_length = static_cast<double>(factors[kYDimension]);
  auto z_length = static_cast<double>(factors[kZDimension]);
  auto basis_of_input = in_config.GetBasis();
  out_config.SetBasis({
      basis_of_input[kXDimension] * x_length,
      basis_of_input[kYDimension] * y_length,
      basis_of_input[kZDimension] * z_length
  });
  out_config.SetScale(1.0);
  auto atom_list_of_input = in_config.GetAtomList();
  int atoms_counter = 0;
  for (int k = 0; k < factors[kZDimension]; ++k) {
    for (int j = 0; j < factors[kYDimension]; ++j) {
      for (int i = 0; i < factors[kXDimension]; ++i) {
        auto x_reference = static_cast<double>(i);
        auto y_reference = static_cast<double>(j);
        auto z_reference = static_cast<double>(k);
        for (const auto &atom : atom_list_of_input) {
          out_config.AppendAtom({
              atoms_counter++, atom.mass_, atom.type_,
              (x_reference + atom.relative_position_[kXDimension]) / x_length,
              (y_reference + atom.relative_position_[kYDimension]) / y_length,
              (z_reference + atom.relative_position_[kZDimension]) / z_length
          });
        }
      }
    }
  }
  return out_config;
}

} // namespace kn
