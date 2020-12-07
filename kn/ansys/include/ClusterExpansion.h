#ifndef KN_KN_ANSYS_INCLUDE_CLUSTEREXPANSION_H_
#define KN_KN_ANSYS_INCLUDE_CLUSTEREXPANSION_H_
#include <utility>
#include <vector>
#include <map>
#include <unordered_set>
#include <unordered_map>

#include <boost/functional/hash.hpp>

#include "Config.h"

namespace ansys::ClusterExpansion {
using Cluster_Map_t = std::vector<std::vector<std::vector<size_t>>>;
bool IsAtomSmallerSymmetrically(const cfg::Atom &lhs, const cfg::Atom &rhs);

// this function rotate and sort the config in a particular way, basically from center to outside
// std::array<std::vector<double>, 2> GetAverageClusterParametersForwardAndBackward(
//     const cfg::Config &config,
//     const std::pair<size_t, size_t> &jump_pair,
//     const std::unordered_map<std::string, double> &type_category_hashmap);

// this function rotates and sorts the config in a particular way, basically from center to outside
// and returns a mapping of it
Cluster_Map_t GetAverageClusterParametersMapping(
    const cfg::Config &config);

std::array<std::vector<double>, 2> GetAverageClusterParametersForwardAndBackwardFromMap(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &jump_pair,
    const std::unordered_map<std::string, double> &type_category_hashmap,
    const Cluster_Map_t &cluster_mapping);

} // namespace ansys::ClusterExpansion

#endif //KN_KN_ANSYS_INCLUDE_CLUSTEREXPANSION_H_
