#ifndef KN_SRC_UTILITY_H_
#define KN_SRC_UTILITY_H_
#include <array>
#include "constants.h"
inline std::array<double, kDimension> CrossProduct(const std::array<double,
                                                                    kDimension> &first,
                                                   const std::array<double,
                                                                    kDimension> &second) {

  double result_x = first[kYDim] * second[kZDim] - first[kZDim] * second[kYDim];
  double result_y = first[kZDim] * second[kXDim] - first[kXDim] * second[kZDim];
  double result_z = first[kXDim] * second[kYDim] - first[kYDim] * second[kXDim];
  return {result_x, result_y, result_z};
}

#endif //KN_SRC_UTILITY_H_
