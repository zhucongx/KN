//
// Created by Zhucong Xi on 2/11/20.
//
#include "ElemInfo.h"

double ElemInfo::findMass(const std::string elem) {
  /*this "it" is indeed a iterator*/
  auto it = std::find(elementList.begin(), elementList.end(), elem);
  int index = std::distance(elementList.begin(), it);
  return massList[index];
}
double ElemInfo::findMass(const int atomicNum) {
  return massList[atomicNum];
}