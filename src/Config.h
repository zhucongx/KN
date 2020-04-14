#ifndef KN_SRC_CONFIG_H_
#define KN_SRC_CONFIG_H_

#include <string>
#include <iostream>
#include <sstream>
#include <array>
#include <vector>

#include "armadillo"

#include "Atom.h"

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
  virtual void UpdateNeighbors(double first_r_cutoff, double second_r_cutoff);
  // update both atoms' relative and absolute positions according to periodic
  // boundary condition
  void WrapPeriodicAtomPosition();
  void ShiftAtomToCentral(const Atom::Rank &id);
  void MoveAbsoluteDistance(const Vector3<double> &distance_vector);
  std::map<std::string, int> CountAllBonds(double r_cutoff);
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
                   const Vector3<int> &factors);
  void GenerateHCP(const double &lattice_constant_a,
                   const double &lattice_constant_c,
                   const std::string &element,
                   const Vector3<int> &factors);
 protected:
  double scale_{};
  // lowx, lowy, lowz, highx, highy, highz, xy xz yz
  // std::array<double, 9> cell;
  // length of three edges
  // Vector3<double> length;
  // This three vectors form a matrix matching the matrix in Config and POSCAR file
  // representing three Bravais lattice vector
  Vector3<double> first_bravais_vector_{};
  Vector3<double> second_bravais_vector_{};
  Vector3<double> third_bravais_vector_{};
  // Three translational Bravais lattice vector
  // DoubleVecfirst_reciprocal_vector_{},
  //     second_reciprocal_vector_{}, third_reciprocal_vector_{};
  int num_atoms_{};
  double energy_{};
  // The index of atom in the vector is always same as of the id of the atom
  std::vector<Atom> atom_list_;
  // indicate if the Config has found Atoms' neighbor list
  bool neighborFound;
  std::map<std::string, std::vector<Atom::Rank>> element_list_set_;
};

}// namespace box

#endif //KN_SRC_CONFIG_H_
