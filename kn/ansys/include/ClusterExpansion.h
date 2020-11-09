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
bool IsAtomSmallerSymmetrically(const cfg::Atom &lhs, const cfg::Atom &rhs);

// this function rotate and sort the config in a particular way, basically from center to outside
std::vector<double> GetAverageClusterParametersForward(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &jump_pair,
    const std::unordered_map<std::string, double> &type_category_hashmap);

std::vector<double> GetAverageClusterParametersBackward(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &jump_pair,
    const std::unordered_map<std::string, double> &type_category_hashmap);

std::vector<std::vector<size_t>>  GetAverageClusterParametersForwardMap(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &jump_pair,
    const std::unordered_map<std::string, double> &type_category_hashmap);

std::vector<std::vector<size_t>>  GetAverageClusterParametersBackwardMap(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &jump_pair,
    const std::unordered_map<std::string, double> &type_category_hashmap);

std::vector<double> GetAverageClusterParametersForwardFromMap(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &jump_pair,
    const std::unordered_map<std::string, double> &type_category_hashmap,
    const std::vector<std::vector<size_t>> &forward_map);

std::vector<double> GetAverageClusterParametersBackwardFromMap(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &jump_pair,
    const std::unordered_map<std::string, double> &type_category_hashmap,
    const std::vector<std::vector<size_t>> &backward_map);

} // namespace ansys::ClusterExpansion

#endif //KN_KN_ANSYS_INCLUDE_CLUSTEREXPANSION_H_
