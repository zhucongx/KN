#ifndef KN_SRC_CONFIG_H_
#define KN_SRC_CONFIG_H_

#include <string>
#include <iostream>
#include <sstream>
#include <array>
#include <vector>

#include "armadillo"

#include "Atom.h"
#include "Box.h"

class Config {
 public:
  /*
  Config.cpp
  */
  void ConvertRelativeToAbsolute();
  void ConvertAbsoluteToRelative();
  Config();
  bool operator<(const Config &rhs) const;
  void Initialize();
  // void ConvertBravisLatticeToReciprocal();
  void Perturb();
  /*
  ConfigIO.cpp
  */
  bool ReadConfig(const std::string &file_name);
  bool ReadPOSCAR(const std::string &file_name);
  void WriteConfig(const std::string &file_name = "config") const;
  // Write Configuration out as POSCAR file. If the show_vacancy_option is
  // true, output will have "X" for visualization. If false, vacancies will be
  // ignored for VASP calculation.
  void WritePOSCAR(const std::string &file_name = "POSCAR",
                   const bool &show_vacancy_option = false) const;
  /*
  ConfigGenerate.cpp
  */
  void GenerateFCC(const double &lattice_constant_a,
                   const std::string &element,
                   const std::array<int, kDimension> &factors);
  void GenerateBCC(const double &lattice_constant_a,
                   const std::string &element,
                   const std::array<int, kDimension> &factors);
  void GenerateHCP(const double &lattice_constant_a,
                   const double &lattice_constant_c,
                   const std::string &element,
                   const std::array<int, kDimension> &factors);
 private:
  int num_atoms_{};
  double energy_{};
  Box box_;
  std::vector<Atom> atom_list_;
  std::vector<int> vacancy_list_;
};

#endif //KN_SRC_CONFIG_H_
