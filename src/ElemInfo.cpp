//
// Created by Zhucong Xi on 2/11/20.
//
#include "ElemInfo.h"

namespace ElemInfo {
double findMass(const std::string elem) {
  /*this "it" is indeed a iterator*/
  auto it = std::find(elementList.begin(), elementList.end(), elem);
  int index = std::distance(elementList.begin(), it);
  return massList[index];
}
double findMass(const int atomicNum) {
  return massList[atomicNum];
}
} // namespace ElemInfo
