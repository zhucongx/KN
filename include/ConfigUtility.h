#ifndef KN_INCLUDE_CONFIGUTILITY_H_
#define KN_INCLUDE_CONFIGUTILITY_H_

#include <vector>
#include <map>
#include "Bond.h"
#include "Config.h"

namespace kn {
class ConfigUtility {
  public:
    static std::map<Bond, int> CountAllBonds(Config &config);
    // add small perturbation to break perfect fcc symmetry this method is about to increase
    // the chance to find lower ground states for VASP software
    static std::vector<std::vector<std::string>> Encode(const Config &config,
                                                        const std::pair<int, int> &jump_pair);

};
} // namespace kn


#endif //KN_INCLUDE_CONFIGUTILITY_H_
