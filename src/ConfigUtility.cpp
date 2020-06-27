#include "ConfigUtility.h"
#include <unordered_set>
namespace kn {

std::map<Bond, int> ConfigUtility::CountAllBonds(Config &config) {
  std::map<Bond, int> bonds_count_map;
  std::string type1, type2;
  auto atom_list = config.GetAtomList();
  for (const auto &atom : atom_list) {
    type1 = atom.type_;
    for (const auto &atom2_id : atom.first_nearest_neighbor_list_) {
      bonds_count_map[Bond{type1, atom_list[atom2_id].type_}]++;
    }
  }
  for (auto &bond_count : bonds_count_map) {
    bond_count.second /= 2;
  }
  return bonds_count_map;
}
} // namespace kn
