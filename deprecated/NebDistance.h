#ifndef KN_INCLUDE_NEBDISTANCE_H_
#define KN_INCLUDE_NEBDISTANCE_H_
#include <vector>
#include <map>
#include <unordered_set>
#include <unordered_map>

#include "Bond.h"
#include "Cfg.hpp"

namespace kn::NebDistance {
void PrintTheDistanceFromTwoPOSCARFiles(const std::string &filename1,
                                        const std::string &filename2);

} // namespace kn::NebDistance

#endif //KN_INCLUDE_NEBDISTANCE_H_
