#include "WarrenCowley.h"
namespace ansys {

WarrenCowley::WarrenCowley(const std::string &cfg_filename)
    : config_(cfg::Config::ReadConfig(cfg_filename, 2)) {
  for (const auto &[element_type, index_vector]: config_.GetElementListMap()) {
    if (element_type == "X")
      continue;
    element_set_.insert(element_type);
  }
}
std::map<std::string, double> WarrenCowley::FindWarrenCowley() const {
  std::map<std::string, double> res;
  for (const auto &type1: element_set_) {
    for (const auto &type2: element_set_) {
      res[std::string(type1).append("-").append(type2)] = 0;
    }
  }
  const auto &element_list_map = config_.GetElementListMap();
  std::map<std::string, double> concentration;
  std::map<std::string, double> count;
  for (const auto &type: element_set_) {
    count[type] = static_cast<double>(element_list_map.at(type).size());
    concentration[type] = static_cast<double>(element_list_map.at(type).size())
        / static_cast<double>(config_.GetNumAtoms());
  }

  for (const auto &atom1: config_.GetAtomList()) {
    const auto type1 = atom1.GetType();
    if (type1 == "X") {
      continue;
    }
    const auto &neighbor_list1 = atom1.GetFirstNearestNeighborsList();
    for (const auto &type2: element_set_) {
      double bond_count = 0;
      for (auto neighbor_index: neighbor_list1) {
        if (config_.GetAtomList()[neighbor_index].GetType() == type2) {
          bond_count += 1;
        }
      }
      double pij = bond_count / 12;
      double aij = (pij - concentration[type2])
          / (static_cast<double>(atom1.GetType() == type2) - concentration[type2]);

      res[std::string(type1).append("-").append(type2)] += aij / count[type1];
    }
  }
  // for (const auto &type1: element_set_) {
  //   for (const auto &type2: element_set_) {
  //     std::cout << std::string(type1).append("-").append(type2) << ' '
  //               << res[std::string(type1).append("-").append(type2)] << '\n';
  //   }
  // }

  return res;

}
} // namespace ansys