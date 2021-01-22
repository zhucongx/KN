#include "ConfigGenerator.h"
#include <chrono>
#include <filesystem>
#include <utility>

#include "InitialConfigGenerator.h"
namespace gen {
InitialConfigGenerator::InitialConfigGenerator(double lattice_const,
                                               const Factor_t &factors,
                                               std::string solvent_element,
                                               std::map<std::string, size_t> element_count_map) :
    lattice_const_(lattice_const),
    factors_(factors),
    solvent_element_(std::move(solvent_element)),
    element_count_map_(std::move(element_count_map)),
    generator_(static_cast<size_t>(std::chrono::system_clock::now().time_since_epoch().count())) {
  std::transform(element_count_map_.begin(),
                 element_count_map_.end(),
                 std::back_inserter(element_list_),
                 [](const auto &pair) { return pair.first; });

  for (const auto&[element, count] : element_count_map_) {
    atom_index_list_.insert(atom_index_list_.end(), count, element);
  }
}
cfg::Config InitialConfigGenerator::EmbedToLarge(double lattice_const,
                                                 const Factor_t &new_factors,
                                                 const Factor_t &old_factors,
                                                 const cfg::Config &small_config,
                                                 std::map<std::string, size_t> element_number_map) {
  size_t new_all = new_factors[0] * new_factors[1] * new_factors[2] * 4;
  size_t old_all = old_factors[0] * old_factors[1] * old_factors[2] * 4;
  std::map<std::string, size_t> old_type_count = cfg::CountAllType(small_config);

  for (const auto &[type, count] : element_number_map) {
    std::cerr << type << element_number_map[type] << '\n';
    element_number_map.at(type) -= old_type_count[type];
    std::cerr << type << old_type_count[type] << '\n';
  }
  std::vector<std::string> atom_type_list{};
  for (const auto&[element, count] : element_number_map) {
    std::cerr << element << count << '\n';

    atom_type_list.insert(atom_type_list.end(), count, element);
  }
  std::mt19937_64
      generator(static_cast<size_t>(std::chrono::system_clock::now().time_since_epoch().count()));
  shuffle(atom_type_list.begin(), atom_type_list.end(), generator);

  cfg::Config config({{{lattice_const * static_cast<double>(new_factors[kXDimension]), 0, 0},
                       {0, lattice_const * static_cast<double>(new_factors[kYDimension]), 0},
                       {0, 0, lattice_const * static_cast<double>(new_factors[kZDimension])}}},
                     small_config.GetAtomList());
  config.ConvertCartesianToRelative();
  size_t atoms_counter = 0;
  std::string element;
  double mass;
  auto x_length = static_cast<double>(new_factors[kXDimension]);
  auto y_length = static_cast<double>(new_factors[kYDimension]);
  auto z_length = static_cast<double>(new_factors[kZDimension]);
  for (size_t k = 0; k < new_factors[kZDimension]; ++k) {
    for (size_t j = 0; j < new_factors[kYDimension]; ++j) {
      for (size_t i = 0; i < new_factors[kXDimension]; ++i) {
        if (k < old_factors[kZDimension] && j < old_factors[kYDimension]
            && i < old_factors[kXDimension])
          continue;
        auto x_reference = static_cast<double>(i);
        auto y_reference = static_cast<double>(j);
        auto z_reference = static_cast<double>(k);
        element = atom_type_list[atoms_counter];
        mass = cfg::FindMass(element);
        config.AppendAtomWithChangingAtomID({
                                                atoms_counter++, mass, element,
                                                x_reference / x_length,
                                                y_reference / y_length,
                                                z_reference / z_length
                                            });
        element = atom_type_list[atoms_counter];
        mass = cfg::FindMass(element);
        config.AppendAtomWithChangingAtomID({
                                                atoms_counter++, mass, element,
                                                (x_reference + 0.5) / x_length,
                                                (y_reference + 0.5) / y_length,
                                                z_reference / z_length
                                            });
        element = atom_type_list[atoms_counter];
        mass = cfg::FindMass(element);
        config.AppendAtomWithChangingAtomID({
                                                atoms_counter++, mass, element,
                                                (x_reference + 0.5) / x_length,
                                                y_reference / y_length,
                                                (z_reference + 0.5) / z_length
                                            });
        element = atom_type_list[atoms_counter];
        mass = cfg::FindMass(element);
        config.AppendAtomWithChangingAtomID({
                                                atoms_counter++, mass, element,
                                                x_reference / x_length,
                                                (y_reference + 0.5) / y_length,
                                                (z_reference + 0.5) / z_length
                                            });
      }
    }
  }
  config.UpdateNeighbors();

  return config;
}
cfg::Config InitialConfigGenerator::ShuffleConfig(const cfg::Config &config) const {
  cfg::Config config_out(config.GetBasis(), config.GetNumAtoms());
  auto atom_index_list_copy(atom_index_list_);
  shuffle(atom_index_list_copy.begin(), atom_index_list_copy.end(), generator_);
  for (auto atom : config.GetAtomList()) {
    atom.SetType(atom_index_list_copy[atom.GetId()]);
    config_out.AppendAtomWithoutChangingAtomID(atom);
  }
  return config_out;
}

static cfg::Config GenerateUnitCell(
    const Matrix_t &basis_matrix,
    const std::vector<std::pair<std::string, Vector_t>> &type_position_list) {
  cfg::Config config(basis_matrix, type_position_list.size());
  size_t atoms_counter = 0;
  for (const auto &[type, relative_position] : type_position_list) {
    config.AppendAtomWithoutChangingAtomID({atoms_counter++, cfg::FindMass(type), type,
                                            relative_position});
  }
  config.ConvertRelativeToCartesian();
  return config;
}

static cfg::Config Duplicate(const cfg::Config &in_config,
                             const Factor_t &factors) {
  auto x_length = static_cast<double>(factors[kXDimension]);
  auto y_length = static_cast<double>(factors[kYDimension]);
  auto z_length = static_cast<double>(factors[kZDimension]);
  auto basis_of_input = in_config.GetBasis();
  cfg::Config out_config({basis_of_input[kXDimension] * x_length,
                          basis_of_input[kYDimension] * y_length,
                          basis_of_input[kZDimension] * z_length},
                         in_config.GetNumAtoms() * factors[kXDimension] * factors[kYDimension]
                             * factors[kZDimension]);
  auto atom_list_of_input = in_config.GetAtomList();
  size_t atoms_counter = 0;
  for (size_t k = 0; k < factors[kZDimension]; ++k) {
    for (size_t j = 0; j < factors[kYDimension]; ++j) {
      for (size_t i = 0; i < factors[kXDimension]; ++i) {
        auto x_reference = static_cast<double>(i);
        auto y_reference = static_cast<double>(j);
        auto z_reference = static_cast<double>(k);
        for (const auto &atom : atom_list_of_input) {
          out_config.AppendAtomWithoutChangingAtomID({atoms_counter++, atom.GetMass(),
                                                      atom.GetType(),
                                                      (x_reference
                                                          + atom.GetRelativePosition()[kXDimension])
                                                          / x_length,
                                                      (y_reference
                                                          + atom.GetRelativePosition()[kYDimension])
                                                          / y_length,
                                                      (z_reference
                                                          + atom.GetRelativePosition()[kZDimension])
                                                          / z_length
                                                     });
        }
      }
    }
  }
  out_config.ConvertRelativeToCartesian();
  return out_config;
}

cfg::Config InitialConfigGenerator::GenerateL10(double lattice_constant_a,
                                                const std::vector<std::string> &element_list,
                                                const Factor_t &factors) {
  Matrix_t basis = {{
                        {lattice_constant_a, 0, 0},
                        {0, lattice_constant_a, 0},
                        {0, 0, lattice_constant_a}
                    }};
  std::vector<std::pair<std::string, Vector_t>> element_position_list;
  element_position_list.emplace_back(element_list[0], Vector_t{0.0, 0.0, 0.0});
  element_position_list.emplace_back(element_list[0], Vector_t{0.5, 0.5, 0.0});
  element_position_list.emplace_back(element_list[1], Vector_t{0.0, 0.5, 0.5});
  element_position_list.emplace_back(element_list[1], Vector_t{0.5, 0.0, 0.5});

  return Duplicate(GenerateUnitCell(basis, element_position_list), factors);
}

cfg::Config InitialConfigGenerator::GenerateL12(double lattice_constant_a,
                                                const std::vector<std::string> &element_list,
                                                const Factor_t &factors) {
  Matrix_t basis = {
      {
          {lattice_constant_a, 0, 0},
          {0, lattice_constant_a, 0},
          {0, 0, lattice_constant_a}
      }
  };
  std::vector<std::pair<std::string, Vector_t>> element_position_list;
  element_position_list.emplace_back(element_list[0], Vector_t{0.0, 0.0, 0.0});
  element_position_list.emplace_back(element_list[1], Vector_t{0.5, 0.5, 0.0});
  element_position_list.emplace_back(element_list[1], Vector_t{0.0, 0.5, 0.5});
  element_position_list.emplace_back(element_list[1], Vector_t{0.5, 0.0, 0.5});

  return Duplicate(GenerateUnitCell(basis, element_position_list), factors);
}

cfg::Config InitialConfigGenerator::GenerateL10star(double lattice_constant_a,
                                                    const std::vector<std::string> &element_list,
                                                    const Factor_t &factors) {
  Matrix_t basis = {
      {
          {lattice_constant_a * 2, 0, 0},
          {0, lattice_constant_a, 0},
          {0, 0, lattice_constant_a}
      }
  };
  std::vector<std::pair<std::string, Vector_t>> element_position_list;
  element_position_list.emplace_back(element_list[0], Vector_t{0.0, 0.0, 0.0});
  element_position_list.emplace_back(element_list[0], Vector_t{0.5, 0.0, 0.0});
  element_position_list.emplace_back(element_list[1], Vector_t{0.25, 0.5, 0.0});
  element_position_list.emplace_back(element_list[0], Vector_t{0.75, 0.5, 0.0});
  element_position_list.emplace_back(element_list[0], Vector_t{0.25, 0.0, 0.5});
  element_position_list.emplace_back(element_list[1], Vector_t{0.75, 0.0, 0.5});
  element_position_list.emplace_back(element_list[1], Vector_t{0.0, 0.5, 0.5});
  element_position_list.emplace_back(element_list[1], Vector_t{0.5, 0.5, 0.5});

  return Duplicate(GenerateUnitCell(basis, element_position_list), factors);
}

cfg::Config InitialConfigGenerator::GenerateL12star(double lattice_constant_a,
                                                    const std::vector<std::string> &element_list,
                                                    const Factor_t &factors) {
  Matrix_t basis = {
      {
          {lattice_constant_a * 4, 0, 0},
          {0, lattice_constant_a, 0},
          {0, 0, lattice_constant_a}
      }
  };
  std::vector<std::pair<std::string, Vector_t>> element_position_list;
  element_position_list.emplace_back(element_list[0], Vector_t{0.0, 0.0, 0.0});
  element_position_list.emplace_back(element_list[1], Vector_t{0.125, 0.5, 0.0});
  element_position_list.emplace_back(element_list[1], Vector_t{0.25, 0.0, 0.0});
  element_position_list.emplace_back(element_list[1], Vector_t{0.375, 0.5, 0.0});
  element_position_list.emplace_back(element_list[1], Vector_t{0.5, 0.0, 0.0});
  element_position_list.emplace_back(element_list[1], Vector_t{0.625, 0.5, 0.0});
  element_position_list.emplace_back(element_list[0], Vector_t{0.75, 0.0, 0.0});
  element_position_list.emplace_back(element_list[1], Vector_t{0.875, 0.5, 0.0});
  element_position_list.emplace_back(element_list[1], Vector_t{0.0, 0.5, 0.5});
  element_position_list.emplace_back(element_list[1], Vector_t{0.125, 0.0, 0.5});
  element_position_list.emplace_back(element_list[0], Vector_t{0.25, 0.5, 0.5});
  element_position_list.emplace_back(element_list[1], Vector_t{0.375, 0.0, 0.5});
  element_position_list.emplace_back(element_list[0], Vector_t{0.5, 0.5, 0.5});
  element_position_list.emplace_back(element_list[1], Vector_t{0.625, 0.0, 0.5});
  element_position_list.emplace_back(element_list[1], Vector_t{0.75, 0.5, 0.5});
  element_position_list.emplace_back(element_list[1], Vector_t{0.875, 0.0, 0.5});

  return Duplicate(GenerateUnitCell(basis, element_position_list), factors);
}

cfg::Config InitialConfigGenerator::GenerateZ1(double lattice_constant_a,
                                               const std::vector<std::string> &element_list,
                                               const Factor_t &factors) {
  Matrix_t basis = {
      {
          {lattice_constant_a, 0, 0},
          {0, lattice_constant_a, 0},
          {0, 0, lattice_constant_a * 2}
      }
  };
  std::vector<std::pair<std::string, Vector_t>> element_position_list;
  element_position_list.emplace_back(element_list[0], Vector_t{0.0, 0.0, 0.0});
  element_position_list.emplace_back(element_list[0], Vector_t{0.5, 0.5, 0.0});
  element_position_list.emplace_back(element_list[1], Vector_t{0.0, 0.5, 0.25});
  element_position_list.emplace_back(element_list[1], Vector_t{0.5, 0.0, 0.25});
  element_position_list.emplace_back(element_list[1], Vector_t{0.0, 0.0, 0.5});
  element_position_list.emplace_back(element_list[1], Vector_t{0.5, 0.5, 0.5});
  element_position_list.emplace_back(element_list[1], Vector_t{0.0, 0.5, 0.75});
  element_position_list.emplace_back(element_list[1], Vector_t{0.5, 0.0, 0.75});

  return Duplicate(GenerateUnitCell(basis, element_position_list), factors);
}

} // namespace cfg
