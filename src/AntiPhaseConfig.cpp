#include "AntiPhaseConfig.h"
namespace box
{
void AntiPhaseConfig::GenerateL10(const double &lattice_constant_a,
                                  const std::vector<std::string> &element_list,
                                  const std::array<int, kDimension> &factors)
{
  Matrix33 bravais_matrix = {{{lattice_constant_a, 0, 0},
                              {0, lattice_constant_a, 0},
                              {0, 0, lattice_constant_a}}};
  std::vector<std::pair<std::string, Vector3>> element_position_list;
  element_position_list.emplace_back(element_list[0], Vector3{0.0, 0.0, 0.0});
  element_position_list.emplace_back(element_list[0], Vector3{0.5, 0.5, 0.0});
  element_position_list.emplace_back(element_list[1], Vector3{0.0, 0.5, 0.5});
  element_position_list.emplace_back(element_list[1], Vector3{0.5, 0.0, 0.5});

  GenerateUnitCell(bravais_matrix, element_position_list);
  Duplicate(factors);
}
void AntiPhaseConfig::GenerateL12(const double &lattice_constant_a,
                                  const std::vector<std::string> &element_list,
                                  const std::array<int, kDimension> &factors)
{
  Matrix33 bravais_matrix = {{{lattice_constant_a, 0, 0},
                              {0, lattice_constant_a, 0},
                              {0, 0, lattice_constant_a}}};
  std::vector<std::pair<std::string, Vector3>> element_position_list;
  element_position_list.emplace_back(element_list[0], Vector3{0.0, 0.0, 0.0});
  element_position_list.emplace_back(element_list[1], Vector3{0.5, 0.5, 0.0});
  element_position_list.emplace_back(element_list[1], Vector3{0.0, 0.5, 0.5});
  element_position_list.emplace_back(element_list[1], Vector3{0.5, 0.0, 0.5});

  GenerateUnitCell(bravais_matrix, element_position_list);
  Duplicate(factors);
}
void AntiPhaseConfig::GenerateL10star(const double &lattice_constant_a,
                                      const std::vector<std::string> &element_list,
                                      const std::array<int, kDimension> &factors)
{
  Matrix33 bravais_matrix = {{{lattice_constant_a * 2, 0, 0},
                              {0, lattice_constant_a, 0},
                              {0, 0, lattice_constant_a}}};
  std::vector<std::pair<std::string, Vector3>> element_position_list;
  element_position_list.emplace_back(element_list[0], Vector3{0.0, 0.0, 0.0});
  element_position_list.emplace_back(element_list[0], Vector3{0.5, 0.0, 0.0});
  element_position_list.emplace_back(element_list[1], Vector3{0.25, 0.5, 0.0});
  element_position_list.emplace_back(element_list[0], Vector3{0.75, 0.5, 0.0});
  element_position_list.emplace_back(element_list[0], Vector3{0.25, 0.0, 0.5});
  element_position_list.emplace_back(element_list[1], Vector3{0.75, 0.0, 0.5});
  element_position_list.emplace_back(element_list[1], Vector3{0.0, 0.5, 0.5});
  element_position_list.emplace_back(element_list[1], Vector3{0.5, 0.5, 0.5});

  GenerateUnitCell(bravais_matrix, element_position_list);
  Duplicate(factors);
}
void AntiPhaseConfig::GenerateL12star(const double &lattice_constant_a,
                                      const std::vector<std::string> &element_list,
                                      const std::array<int, kDimension> &factors)
{
  Matrix33 bravais_matrix = {{{lattice_constant_a * 4, 0, 0},
                              {0, lattice_constant_a, 0},
                              {0, 0, lattice_constant_a}
                             }};
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
  GenerateUnitCell(bravais_matrix, element_position_list);
  Duplicate(factors);
}

void AntiPhaseConfig::GenerateZ1(const double &lattice_constant_a,
                                 const std::vector<std::string> &element_list,
                                 const std::array<int, kDimension> &factors)
{
  Matrix33 bravais_matrix = {{{lattice_constant_a, 0, 0},
                              {0, lattice_constant_a, 0},
                              {0, 0, lattice_constant_a * 2}}};
  std::vector<std::pair<std::string, Vector3>> element_position_list;
  element_position_list.emplace_back(element_list[0], Vector3{0.0, 0.0, 0.0});
  element_position_list.emplace_back(element_list[0], Vector3{0.5, 0.5, 0.0});
  element_position_list.emplace_back(element_list[1], Vector3{0.0, 0.5, 0.25});
  element_position_list.emplace_back(element_list[1], Vector3{0.5, 0.0, 0.25});
  element_position_list.emplace_back(element_list[1], Vector3{0.0, 0.0, 0.5});
  element_position_list.emplace_back(element_list[1], Vector3{0.5, 0.5, 0.5});
  element_position_list.emplace_back(element_list[1], Vector3{0.0, 0.5, 0.75});
  element_position_list.emplace_back(element_list[1], Vector3{0.5, 0.0, 0.75});

  GenerateUnitCell(bravais_matrix, element_position_list);
  Duplicate(factors);
}

}// namespace box