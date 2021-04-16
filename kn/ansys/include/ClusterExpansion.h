#ifndef KN_KN_ANSYS_INCLUDE_CLUSTEREXPANSION_H_
#define KN_KN_ANSYS_INCLUDE_CLUSTEREXPANSION_H_
#include <utility>
#include <vector>
#include <map>
#include <unordered_set>
#include <unordered_map>

#include "Config.h"

namespace ansys::ClusterExpansion {
std::unordered_map<std::string, std::vector<double> > GetOneHotEncodeHashmap(
    const std::set<std::string> &type_set);

// this function rotates and sorts the config in a particular way, basically from center to outside
// and returns a mapping of it
std::vector<std::vector<std::vector<size_t> > > GetAverageClusterParametersMapping(const cfg::Config &config);

std::array<std::vector<std::string>, 2> GetForwardAndBackwardEncode(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &atom_id_jump_pair);

std::vector<double> GetOneHotParametersFromMap(
    const std::vector<std::string> &encode,
    const std::unordered_map<std::string, std::vector<double> > &one_hot_encode_hashmap,
    size_t num_of_elements,
    const std::vector<std::vector<std::vector<size_t> > > &cluster_mapping);

} // namespace ansys::ClusterExpansion

#endif //KN_KN_ANSYS_INCLUDE_CLUSTEREXPANSION_H_
