#ifndef KN_SRC_CONFIG_H_
#define KN_SRC_CONFIG_H_

#include <string>
#include <iostream>
#include <sstream>
#include <array>
#include <vector>

#include "armadillo"

#include "Atom.h"
#include "Cell.h"
#include "Utility.h"

namespace box {

class Config {
 public:
  /*
  Config.cpp
  */
  Config();
  bool operator<(const Config &rhs) const;
  void Perturb();
  virtual void UpdateNeighbors(double firrst_r_cutoff, double second_r_cutoff);
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

 protected:
  /*
  Config.cpp
  */
  void Initialize();
  void ConvertRelativeToAbsolute();
  void ConvertAbsoluteToRelative();
  [[nodiscard]] Double3 GetRelativeDistanceVector(int first, int second) const;
  int num_atoms_{};
  double energy_{};
  Cell box_;
  std::vector<Atom> atom_list_;
  std::vector<int> vacancy_list_;
};

}// namespace box

#endif //KN_SRC_CONFIG_H_
