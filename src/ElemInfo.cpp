
#include "ElemInfo.h"

namespace elem_info {
double FindMass(const std::string &elem) {
  /*this "it" is indeed a iterator*/
  auto it = std::find(element_list.begin(), element_list.end(), elem);
  return mass_list[std::distance(element_list.begin(), it)];
}
double FindMass(const int &atomicNum) {
  return mass_list[atomicNum];
}
} // namespace elem_info
