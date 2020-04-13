#ifndef KN_SRC_HCPCONFIG_H_
#define KN_SRC_HCPCONFIG_H_

#include "Config.h"

namespace box {

class HCPConfig : public Config {
 public:
  void GenerateHCP(const double &lattice_constant_a,
                   const double &lattice_constant_c,
                   const std::string &element,
                   const Int3 &factors);
};

}// namespace box

#endif //KN_SRC_HCPCONFIG_H_
