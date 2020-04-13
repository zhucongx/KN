#ifndef KN_SRC_BCCCONFIG_H_
#define KN_SRC_BCCCONFIG_H_

#include "Config.h"

namespace box {

class BCCConfig : public Config {
 public:
  void GenerateBCC(const double &lattice_constant_a,
                   const std::string &element,
                   const Int3 &factors);
  void UpdateNeighbors(double first_r_cutoff, double second_r_cutoff) override;
};

}// namespace box

#endif //KN_SRC_BCCCONFIG_H_
