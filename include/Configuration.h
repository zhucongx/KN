//
// Created by Zhucong Xi on 1/31/20.
//

#ifndef _CONFIGURATION_H_
#define _CONFIGURATION_H_

#include <string>
#include <fstream>
#include <iostream>
#include <array>
#include <vector>
#include "Atom.h"

using std::array;
using std::vector;
using std::ifstream;

enum { XX = 0, YY = 1, ZZ = 2, XY = 3, YZ = 4, ZX = 5 };
enum { X = 0, Y = 1, Z = 2 };

class Configuration {
 private:
  int numAtoms;
  int numTypes;
  double energy;
  // lowx, lowy, lowz, highx, highy, highz, xy xz yz
  array<double, 9> cell;
  // length of three edges
  array<double, 3> length;
  // bvx, bvy, bvz form a matrix matching H0 in .cfg file
  array<double, 3> bvx, tvx, bvy, tvy, bvz, tvz;
  vector<Atom> atoms;
  vector<int> vacList;
 public:
  Configuration();
  ~Configuration();
  bool operator<(const Configuration &rhs) const;
  bool operator>(const Configuration &rhs) const;
  bool operator<=(const Configuration &rhs) const;
  bool operator>=(const Configuration &rhs) const;

  void readLammpsData(const string &fileName);
  bool readConfig(const string &fileName);
  void readPOSCAR(const string &fileName);
};

#endif //KN_SRC_CONFIGURATION_H_
