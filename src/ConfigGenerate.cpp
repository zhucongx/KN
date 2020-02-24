//
// Created by Zhucong Xi on 2/11/20.
//

#include "Config.h"

void Config::generateFCC(const double &latticeConstant,
                         const std::string &elm,
                         const std::vector<int> &factors) {
  clear();
  double ms = ElemInfo::findMass(elm);
  bvx[X] = latticeConstant * factors[X];
  bvy[Y] = latticeConstant * factors[Y];
  bvz[Z] = latticeConstant * factors[Z];
  int c = 0;
  auto na = static_cast<double>(factors[X]);
  auto nb = static_cast<double>(factors[Y]);
  auto nc = static_cast<double>(factors[Z]);
  for (int k = 0; k < factors[Z]; ++k) {
    for (int j = 0; j < factors[Y]; ++j) {
      for (int i = 0; i < factors[X]; ++i) {
        auto x = static_cast<double>(i);
        auto y = static_cast<double>(j);
        auto z = static_cast<double>(k);
        atoms.emplace_back(c++, ms, elm, x / na, y / nb, z / nc);
        atoms.emplace_back(c++, ms, elm, (x + .5) / na, (y + .5) / nb, z / nc);
        atoms.emplace_back(c++, ms, elm, (x + .5) / na, y / nb, (z + .5) / nc);
        atoms.emplace_back(c++, ms, elm, (x + .5) / na, y / nb, (z + .5) / nc);
      }
    }
  }
  numAtoms = c;
}

void Config::embedCluster(const std::pair<std::string, std::string> &Elems,
                          const FCCEmbededCluster::occupInfo_256 &o256,
                          const int &i) {
  for (int j = 0; j < o256.mapping[i].size(); ++j) {
    if (o256.mapping[i][j] == 1) {
      atoms[j].setType(Elems.first) ;
    } else if (o256.mapping[i][j] == 2) {
      atoms[j].setType(Elems.second);
    }
  }
}