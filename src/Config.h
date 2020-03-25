#ifndef KN_SRC_CONFIG_H_
#define KN_SRC_CONFIG_H_

#include <string>
#include <iostream>
#include <sstream>
#include <array>
#include <vector>

#include "Atom.h"

class Config {
 public:
  // Config.cpp
  void ConvertRelativeToAbsolute();
  void ConvertAbsoluteToRelative();
  Config();
  virtual ~Config();
  bool operator<(const Config &rhs) const;
  void clear();
  bool ReadConfig(const std::string &file_name);
  bool ReadPOSCAR(const std::string &file_name);
  void WriteConfig(const std::string &file_name = "config") const;
  // Write Configuration out as POSCAR file. If the show_vacancy_option is
  // true, output will have "X" for visualization. If false, vacancies will be
  // ignored for VASP calculation.
  void WritePOSCAR(const std::string &file_name = "POSCAR",
                   const bool &show_vacancy_option = false) const;

  void generateFCC(const double &latticeConstant, const std::string &elm,
                   const std::vector<int> &factors);
 private:
  int num_atoms_{};
  double energy_{};
  // lowx, lowy, lowz, highx, highy, highz, xy xz yz
  // std::array<double, 9> cell;
  // length of three edges
  // std::array<double, 3> length;
  // bvx, bvy, bvz form a matrix matching the matrix in Config and POSCAR file
  // representing three Bravais lattice vector
  std::array<double, kDimension> x_bravais_vector_{}, y_bravais_vector_{},
      z_bravais_vector_{};
  // Three translational Bravais lattice vector
  // std::array<double, 3> tvx, tvy, tvz;
  std::vector<Atom> atom_list_;
  std::vector<int> vacancy_list_;
};

#endif //KN_SRC_CONFIG_H_
