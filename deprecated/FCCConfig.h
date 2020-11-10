#ifndef KN_INCLUDE_FCCCONFIG_H_
#define KN_INCLUDE_FCCCONFIG_H_

#include "Cfg.hpp"

namespace kn {

class FCCConfig : public Config {
  public:
    void UpdateNeighbors(double first_r_cutoff, double second_r_cutoff) override;
};

} // namespace kn

#endif //KN_INCLUDE_FCCCONFIG_H_
