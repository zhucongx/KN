#include "ConfigGenerator.h"
#include <chrono>
#include <filesystem>
#include <utility>
#include "ConfigIO.h"
namespace kn {
ConfigGenerator::ConfigGenerator(double lattice_const,
                                 const ConfigGenerator::Factor_t &factors,
                                 std::string solvent_element,
                                 std::map<std::string, int> element_count_map,
                                 std::string pot_folder_path) :
    lattice_const_(lattice_const),
    factors_(factors),
    solvent_element_(std::move(solvent_element)),
    element_count_map_(std::move(element_count_map)),
    pot_folder_path_(std::move(pot_folder_path)),
    generator_(std::chrono::system_clock::now().time_since_epoch().count()) {
  std::transform(element_count_map_.begin(),
                 element_count_map_.end(),
                 std::back_inserter(element_list_),
                 [](const auto &pair) { return pair.first; });

  for (const auto&[element, count] : element_count_map_) {
    atom_index_list_.insert(atom_index_list_.end(), count, element);
  }
}
Config ConfigGenerator::ShuffleConfig(const Config &config) const {
  Config config_out(config.GetBasis(), config.GetNumAtoms());

  auto atom_index_list_copy(atom_index_list_);
  shuffle(atom_index_list_copy.begin(), atom_index_list_copy.end(), generator_);
  for (auto atom : config.GetAtomList()) {
    atom.type_ = atom_index_list_copy[atom.id_];
    config_out.AppendAtom(atom);
  }
  return config_out;
}
void ConfigGenerator::CreateRandom(int num_configs) {
  auto config_start = ShuffleConfig(GenerateFCC(lattice_const_, solvent_element_, factors_));

  for (int i = 0; i < num_configs; i++) {
    std::string start_path("config" + std::to_string(i) + "/s");
    std::filesystem::create_directories(start_path);
    ConfigIO::WriteConfig(config_start, start_path + "/start.cfg", false);
    config_start.Perturb(generator_);
    ConfigIO::WritePOSCAR(config_start, "config" + std::to_string(i) + "/s/POSCAR", false);
    PrepareVASPFiles(start_path);
  }
  // const Config config_start = config_start
}
void ConfigGenerator::CreateSpecific() {
  Config config_start = GenerateFCC(lattice_const_, solvent_element_, factors_);

}

Config GenerateUnitCell(
    const Matrix33 &basis_matrix,
    const std::vector<std::pair<std::string, Vector3>> &type_position_list) {
  Config config(basis_matrix, type_position_list.size());
  int atoms_counter = 0;
  for (const auto &[type, relative_position] : type_position_list) {
    config.AppendAtom({atoms_counter++, elem_info::FindMass(type), type, relative_position});
  }
  config.ConvertRelativeToCartesian();
  return config;
}

Config Duplicate(const Config &in_config,
                 const std::array<int, kDimension> &factors) {
  auto x_length = static_cast<double>(factors[kXDimension]);
  auto y_length = static_cast<double>(factors[kYDimension]);
  auto z_length = static_cast<double>(factors[kZDimension]);
  auto basis_of_input = in_config.GetBasis();
  Config out_config({basis_of_input[kXDimension] * x_length,
                     basis_of_input[kYDimension] * y_length,
                     basis_of_input[kZDimension] * z_length},
                    in_config.GetNumAtoms() * factors[kXDimension] * factors[kYDimension]
                        * factors[kZDimension]);
  auto atom_list_of_input = in_config.GetAtomList();
  int atoms_counter = 0;
  for (int k = 0; k < factors[kZDimension]; ++k) {
    for (int j = 0; j < factors[kYDimension]; ++j) {
      for (int i = 0; i < factors[kXDimension]; ++i) {
        auto x_reference = static_cast<double>(i);
        auto y_reference = static_cast<double>(j);
        auto z_reference = static_cast<double>(k);
        for (const auto &atom : atom_list_of_input) {
          out_config.AppendAtom({atoms_counter++, atom.mass_, atom.type_,
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

Config ConfigGenerator::GenerateFCC(double lattice_constant_a,
                                    const std::string &element,
                                    const std::array<int, kDimension> &factors) {

  double mass = elem_info::FindMass(element);
  Config config({{{lattice_constant_a * factors[kXDimension], 0, 0},
                  {0, lattice_constant_a * factors[kYDimension], 0},
                  {0, 0, lattice_constant_a * factors[kZDimension]}}},
                4 * factors[kXDimension] * factors[kYDimension] * factors[kZDimension]);
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

// Config ConfigGenerator::GenerateBCC(double lattice_constant_a,
//                                     const std::string &element,
//                                     const std::array<int, kDimension> &factors) {
//   double mass = elem_info::FindMass(element);
//   Config config({{{lattice_constant_a * factors[kXDimension], 0, 0},
//                   {0, lattice_constant_a * factors[kYDimension], 0},
//                   {0, 0, lattice_constant_a * factors[kZDimension]}}},
//                 2 * factors[kXDimension] * factors[kYDimension] * factors[kZDimension]);
//   int atoms_counter = 0;
//   auto x_length = static_cast<double>(factors[kXDimension]);
//   auto y_length = static_cast<double>(factors[kYDimension]);
//   auto z_length = static_cast<double>(factors[kZDimension]);
//   for (int k = 0; k < factors[kZDimension]; ++k) {
//     for (int j = 0; j < factors[kYDimension]; ++j) {
//       for (int i = 0; i < factors[kXDimension]; ++i) {
//         auto x_reference = static_cast<double>(i);
//         auto y_reference = static_cast<double>(j);
//         auto z_reference = static_cast<double>(k);
//         config.AppendAtom({
//                               atoms_counter++, mass, element,
//                               x_reference / x_length,
//                               y_reference / y_length,
//                               z_reference / z_length
//                           });
//         config.AppendAtom({
//                               atoms_counter++, mass, element,
//                               (x_reference + 0.5) / x_length,
//                               (y_reference + 0.5) / y_length,
//                               (z_reference + 0.5) / z_length
//                           });
//       }
//     }
//   }
//   config.ConvertRelativeToCartesian();
//   return config;
// }
//
Config ConfigGenerator::GenerateHCP(double lattice_constant_a,
                                    double lattice_constant_c,
                                    const std::string &element,
                                    const std::array<int, kDimension> &factors) {

  double mass = elem_info::FindMass(element);
  Config config({{{lattice_constant_a * factors[kXDimension], 0, 0},
                  {-0.5 * lattice_constant_a * factors[kYDimension],
                   0.5 * sqrt(3) * lattice_constant_a * factors[kYDimension], 0},
                  {0, 0, lattice_constant_c * factors[kZDimension]}}},
                2 * factors[kXDimension] * factors[kYDimension] * factors[kZDimension]);
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
  Matrix33 basis = {{
                        {lattice_constant_a, 0, 0},
                        {0, lattice_constant_a, 0},
                        {0, 0, lattice_constant_a}
                    }};
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

} // namespace kn
