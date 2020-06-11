#ifndef KN_INCLUDE_CONFIGUTILITY_H_
#define KN_INCLUDE_CONFIGUTILITY_H_

#include <random>
#include "Bond.h"
#include "Config.h"

namespace kn {
class ConfigUtility {
  public:
    static std::map<Bond, int> CountAllBonds(Config &config, double r_cutoff);
    // add small perturbation to break perfect fcc symmetry this method is about to increase
    // the chance to find lower ground states for VASP software
    static void Perturb(Config &config, std::mt19937_64 &generator);
};
} // namespace kn


#endif //KN_INCLUDE_CONFIGUTILITY_H_
