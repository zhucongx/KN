//
// Created by Zhucong Xi on 1/31/20.
//

#ifndef _CONFIG_H_
#define _CONFIG_H_

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <numeric>
#include <array>
#include <vector>
#include "armadillo"
#include "FCCEmbededCluster.h"
#include "Atom.h"


class Config {
 private:
  int numAtoms;
  double energy;
  // lowx, lowy, lowz, highx, highy, highz, xy xz yz
  // std::array<double, 9> cell;
  // length of three edges
  // std::array<double, 3> length;
  // bvx, bvy, bvz form a matrix matching the matrix in Config and POSCAR file
  // representing three Bravais lattice vector
  std::array<double, 3> bvx, bvy, bvz;
  // Three translational Bravais lattice vector
  // std::array<double, 3> tvx, tvy, tvz;
  std::vector<Atom> atoms;
  // std::vector<int> vacList;
  void cnvPrl2Pst();
  void cnvPst2Prl();
 public:
  Config();
  virtual ~Config();
  bool operator<(const Config &rhs) const;
  void clear();
  bool readConfig(const std::string &fileName);
  bool readPOSCAR(const std::string &fileName);
  void writeConfig(const std::string &fileName = "config") const;
  // Write Configuration out as POSCAR file. If the vacOption is true, output
  // will have "X" for visualization. If false, vacancies will be ignored for
  // VASP calculation.
  void writePOSCAR(const std::string &fileName = "POSCAR",
                   const bool &vacOption = false) const;
  /**
   ConfigGenerate.cpp
   **/
  void generateFCC(const double &latticeConstant, const std::string &elm,
                   const std::vector<int> &factors);
  void embedCluster(const std::pair<std::string, std::string> &Elems,
                    const FCCEmbededCluster::occupInfo_256 &o256, const int &i);
};

#endif //_CONFIG_H_
