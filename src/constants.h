#ifndef KN_SRC_CONSTANTS_H_
#define KN_SRC_CONSTANTS_H_

// 3D simulation
const int kDimension = 3;
// three indexes
const int kXDim = 0, kYDim = 1, kZDim = 2;

struct Double3 {
  double x, y, z;
};
struct Int3 {
  int x, y, z;
};

namespace FCC_const {
const double kRCutoff = 3.5;
// For FCC the first nearest neighbor is 12
const int kNumFirstNearestNeighbors = 12;
// For FCC the first and second nearest neighbor is 18
const int kNumSecondNearNeighbors = 18;

const int kFirstNearestNeighborsReserveSize = 25;
const int kSecondNearestNeighborsReserveSize = 37;
}// namespace FCC_const

#endif //KN_SRC_CONSTANTS_H_

