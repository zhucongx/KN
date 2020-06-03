#ifndef KN_INCLUDE_CONFIGUTILITY_H_
#define KN_INCLUDE_CONFIGUTILITY_H_

#include "Config.h"
namespace kn {
class ConfigUtility {
 public:
  static std::map<Bond, int> CountAllBonds(Config &config, double r_cutoff);

};
}// namespace kn


#endif //KN_INCLUDE_CONFIGUTILITY_H_
