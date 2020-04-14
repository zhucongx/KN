#ifndef KN_SRC_CONFIG_H_
#define KN_SRC_CONFIG_H_

#include <string>
#include <iostream>
#include <sstream>
#include <array>
#include <vector>

#include "armadillo"

#include "Atom.h"
#include "Utility.h"

namespace box {

class Config {
  ///To do finsih GenerateUnitCell
 public:
  /*
  Config.cpp
  */
  Config();
  bool operator<(const Config &rhs) const;
  void Initialize();
  [[nodiscard]] bool IsCubic() const;
  void ConvertRelativeToAbsolute();
  void ConvertAbsoluteToRelative();
  void Perturb();
  virtual void UpdateNeighbors(double firrst_r_cutoff, double second_r_cutoff);
  virtual void ShiftAtomToCentral(const Rank &id);
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
  void GenerateUnitCell(std::vector<std::pair<std::string, Double3>> type_position);
  void GenerateBCC(const double &lattice_constant_a,
                   const std::string &element,
                   const Int3 &factors);
  void GenerateFCC(const double &lattice_constant_a,
                   const std::string &element,
                   const Int3 &factors);
  void GenerateHCP(const double &lattice_constant_a,
                   const double &lattice_constant_c,
                   const std::string &element,
                   const Int3 &factors);
 protected:
  double scale_{};
  // lowx, lowy, lowz, highx, highy, highz, xy xz yz
  // std::array<double, 9> cell;
  // length of three edges
  // Double3 length;
  // This three vectors form a matrix matching the matrix in Config and POSCAR file
  // representing three Bravais lattice vector
  Double3 first_bravais_vector_{};
  Double3 second_bravais_vector_{};
  Double3 third_bravais_vector_{};
  // Three translational Bravais lattice vector
  // DoubleVecfirst_reciprocal_vector_{},
  //     second_reciprocal_vector_{}, third_reciprocal_vector_{};
  int num_atoms_{};
  double energy_{};
  // The index of atom in the vector is always same as of the id of the atom
  std::vector<Atom> atom_list_;
  std::vector<Rank> vacancy_list_;
};

}// namespace box

#endif //KN_SRC_CONFIG_H_
