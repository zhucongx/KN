#ifndef KN_INCLUDE_CLUSTEREXPANSION_H_
#define KN_INCLUDE_CLUSTEREXPANSION_H_

#include "Bond.h"
#include "Config.h"

namespace kn::ClusterExpansion {
// initial and final state
// this function rotate and sort the config in a particular way, basically from center to outside
Config GetStableStateConfig(const Config & config);

// transition state
Config GetTransitionStateConfig(const Config & config);

} // namespace kn::ClusterExpansion

#endif //KN_INCLUDE_CLUSTEREXPANSION_H_
