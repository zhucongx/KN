#include "BondCounting.h"
#include <unordered_map>

namespace ansys::BondCounting {
std::unordered_map<cfg::Bond, size_t, boost::hash<cfg::Bond> > InitializeHashMap(
    const std::set<std::string> &type_set) {
  std::unordered_map<cfg::Bond, size_t, boost::hash<cfg::Bond> > bond_count_hashmap;
  for (size_t label = 1; label <= 7; ++label) {
    for (const auto &element1 : type_set) {
      for (const auto &element2 : type_set) {
        bond_count_hashmap[cfg::Bond(label, element1, element2)] = 0;
      }
    }
  }
  return bond_count_hashmap;
}
std::vector<double> GetBondChange(const cfg::Config &config,
                                  const std::pair<size_t, size_t> &atom_id_jump_pair,
                                  std::unordered_map<cfg::Bond, size_t,
                                                     boost::hash<cfg::Bond> > initialized_hashmap) {
  const auto &atom_list = config.GetAtomList();
  const std::string type1 = atom_list[atom_id_jump_pair.second].GetType();
  // plus new bonds
  for (const auto &atom2_id : atom_list[atom_id_jump_pair.first].GetFirstNearestNeighborsList()) {
    if (atom2_id == atom_id_jump_pair.second)
      continue;
    initialized_hashmap[cfg::Bond{1, type1, atom_list[atom2_id].GetType()}]++;
  }
  for (const auto &atom2_id : atom_list[atom_id_jump_pair.first].GetSecondNearestNeighborsList())
    initialized_hashmap[cfg::Bond{2, type1, atom_list[atom2_id].GetType()}]++;
  for (const auto &atom2_id : atom_list[atom_id_jump_pair.first].GetThirdNearestNeighborsList())
    initialized_hashmap[cfg::Bond{3, type1, atom_list[atom2_id].GetType()}]++;
  for (const auto &atom2_id : atom_list[atom_id_jump_pair.first].GetFourthNearestNeighborsList())
    initialized_hashmap[cfg::Bond{4, type1, atom_list[atom2_id].GetType()}]++;
  for (const auto &atom2_id : atom_list[atom_id_jump_pair.first].GetFifthNearestNeighborsList())
    initialized_hashmap[cfg::Bond{5, type1, atom_list[atom2_id].GetType()}]++;
  for (const auto &atom2_id : atom_list[atom_id_jump_pair.first].GetSixthNearestNeighborsList())
    initialized_hashmap[cfg::Bond{6, type1, atom_list[atom2_id].GetType()}]++;
  for (const auto &atom2_id : atom_list[atom_id_jump_pair.first].GetSeventhNearestNeighborsList())
    initialized_hashmap[cfg::Bond{7, type1, atom_list[atom2_id].GetType()}]++;
  // minus old bonds
  for (const auto &atom2_id : atom_list[atom_id_jump_pair.second].GetFirstNearestNeighborsList()) {
    if (atom2_id == atom_id_jump_pair.first)
      continue;
    initialized_hashmap[cfg::Bond{1, type1, atom_list[atom2_id].GetType()}]--;
  }
  for (const auto &atom2_id : atom_list[atom_id_jump_pair.second].GetSecondNearestNeighborsList())
    initialized_hashmap[cfg::Bond{2, type1, atom_list[atom2_id].GetType()}]--;
  for (const auto &atom2_id : atom_list[atom_id_jump_pair.second].GetThirdNearestNeighborsList())
    initialized_hashmap[cfg::Bond{3, type1, atom_list[atom2_id].GetType()}]--;
  for (const auto &atom2_id : atom_list[atom_id_jump_pair.second].GetFourthNearestNeighborsList())
    initialized_hashmap[cfg::Bond{4, type1, atom_list[atom2_id].GetType()}]--;
  for (const auto &atom2_id : atom_list[atom_id_jump_pair.second].GetFifthNearestNeighborsList())
    initialized_hashmap[cfg::Bond{5, type1, atom_list[atom2_id].GetType()}]--;
  for (const auto &atom2_id : atom_list[atom_id_jump_pair.second].GetSixthNearestNeighborsList())
    initialized_hashmap[cfg::Bond{6, type1, atom_list[atom2_id].GetType()}]--;
  for (const auto &atom2_id : atom_list[atom_id_jump_pair.second].GetSeventhNearestNeighborsList())
    initialized_hashmap[cfg::Bond{7, type1, atom_list[atom2_id].GetType()}]--;

  std::map<cfg::Bond, size_t> ordered(initialized_hashmap.begin(), initialized_hashmap.end());
  std::vector<double> res;
  res.reserve(ordered.size());
  for (const auto &bond_count : ordered) {
    res.push_back(bond_count.second);
  }

  return res;
}
} // namespace ansys::BondCounting