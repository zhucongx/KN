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
};
} // namespace kn


#endif //KN_INCLUDE_CONFIGUTILITY_H_
