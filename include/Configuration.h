//
// Created by Zhucong Xi on 1/31/20.
//

#ifndef _CONFIGURATION_H_
#define _CONFIGURATION_H_


#include <string>
#include <fstream>
#include <array>
#include <vector>
#include "Atom.h"

using std::array;
using std::vector;

class Configuration {
 private:
  int numAtoms;
  int numTypes;
  double energy;
  // lowx, lowy, lowz, highx, highy, highz, xy xz yz
  array<double, 9> cell;
  // length of three edges
  array<double, 3> length;
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

  void readLammpsData(const string & fileName);
  void readConfig(const string & fileName);
  void readPOSCAR(const string & fileName);

};

#endif //KN_SRC_CONFIGURATION_H_
