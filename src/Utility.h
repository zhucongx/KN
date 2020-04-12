#ifndef KN_SRC_UTILITY_H_
#define KN_SRC_UTILITY_H_

#include "constants.h"
inline Double3 CrossProduct(const Double3 &first, const Double3 &second) {

  double result_x = first[kYDim] * second[kZDim] - first[kZDim] * second[kYDim];
  double result_y = first[kZDim] * second[kXDim] - first[kXDim] * second[kZDim];
  double result_z = first[kXDim] * second[kYDim] - first[kYDim] * second[kXDim];
  return {result_x, result_y, result_z};
}

inline double DotProduct(const Double3 &first, const Double3 &second) {
  return first[kXDim] * second[kXDim] + first[kYDim] * second[kYDim]
      + first[kZDim] * second[kZDim];
}
inline Double3 LinearTransform(const Double3 &left,
                               const Double3 &right1,
                               const Double3 &right2,
                               const Double3 &right3) {
  Double3 result;
  for (const auto &i : {kXDim, kYDim, kZDim}) {
    result[i] = left[kXDim] * right1[i] + left[kYDim] * right2[i]
        + left[kZDim] * right3[i];
  }
  return result;
}
#endif //KN_SRC_UTILITY_H_
