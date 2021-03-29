#ifndef KN_KN_ANSYS_INCLUDE_BONDCOUNTING_H_
#define KN_KN_ANSYS_INCLUDE_BONDCOUNTING_H_
#include "Bond.hpp"
#include "Config.h"
namespace ansys::BondCounting {
std::unordered_map<cfg::Bond, int, boost::hash<cfg::Bond>> InitializeHashMap(
    const std::set<std::string> &type_set);
std::vector<double> GetBondChange(const cfg::Config &config,
                                  const std::pair<size_t, size_t> &atom_id_jump_pair,
                                  std::unordered_map<cfg::Bond, int,
                                                     boost::hash<cfg::Bond>> initialized_hashmap);
} // namespace ansys::BondCounting

#endif //KN_KN_ANSYS_INCLUDE_BONDCOUNTING_H_
