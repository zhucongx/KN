//
// Created by Zhucong Xi on 1/31/20.
//

#ifndef _CONFIGURATION_H_
#define _CONFIGURATION_H_

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <numeric>
#include <array>
#include <vector>
#include "armadillo"
#include "Atom.h"

class Configuration {
 private:
  int numAtoms;
  int numTypes;
  double energy;
  // lowx, lowy, lowz, highx, highy, highz, xy xz yz
  std::array<double, 9> cell{};
  // length of three edges
  std::array<double, 3> length{};
  // bvx, bvy, bvz form a matrix matching the matrix in Config and POSCAR file
  // representing three Bravais lattice vector
  std::array<double, 3> bvx{}, bvy{}, bvz{}, tvx{}, tvy{}, tvz{};
  std::vector<Atom> atoms;
  std::vector<int> vacList;
 public:
  Configuration();
  virtual ~Configuration();
  bool operator<(const Configuration &rhs) const;
  bool operator>(const Configuration &rhs) const;
  bool operator<=(const Configuration &rhs) const;
  bool operator>=(const Configuration &rhs) const;

  void cnvPrl2Pst();
  void cnvPst2Prl();

  bool readLammpsData(const std::string &fileName);
  bool readConfig(const std::string &fileName);
  bool readPOSCAR(const std::string &fileName);

  void writeConfig(const std::string &fileName) const;
  // Write Configuration out as POSCAR file. If the vacOption is true, output
  // will have "X" for visualization. If false, vacancies will be ignored for
  // VASP calculation.
  void writePOSCAR(const std::string &fileName, const bool &vacOption) const;
};

#endif //KN_SRC_CONFIGURATION_H_
