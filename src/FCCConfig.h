#ifndef KN_SRC_FCCCONFIG_H_
#define KN_SRC_FCCCONFIG_H_

#include "Config.h"

namespace box {

class FCCConfig : public Config {
 public:
  void UpdateNeighbors(double first_r_cutoff, double second_r_cutoff);
};

}// namespace box

#endif //KN_SRC_FCCCONFIG_H_
