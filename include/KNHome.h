//
// Created by Zhucong Xi on 2/11/20.
//

#ifndef _KNHOME_H_
#define _KNHOME_H_

#include <string>
#include <vector>
#include "FCCEmbededCluster.h"
#include "Config.h"

class KNHome {
 private:
  Config cfg;

 public:
  int generateOrderedSingle(const int &i, int index,
                            const std::vector<int> &dupFactors,
                            const double &latticeConstant,
                            const std::string &POT,
                            const FCCEmbededCluster::occupInfo_256 &o256);


};

#endif //_KNHOME_H_
