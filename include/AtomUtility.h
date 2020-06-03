#ifndef KN_INCLUDE_ATOMUTILITY_H_
#define KN_INCLUDE_ATOMUTILITY_H_

#include "Config.h"
namespace kn {
class AtomUtility {
 public:
  static Vector3 GetRelativeDistanceVector(const Atom &first, const Atom &second);
};
}// namespace kn


#endif //KN_INCLUDE_ATOMUTILITY_H_
