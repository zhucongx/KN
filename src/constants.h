#ifndef KN_SRC_CONSTANTS_H_
#define KN_SRC_CONSTANTS_H_

// 3D simulation
const int kDimension = 3;
// three indexes
const int kXDim = 0, kYDim = 1, kZDim = 2;

typedef std::array<double,kDimension> Double3;
typedef std::array<int,kDimension> Int3;


namespace FCC_const {
const double kRCutoff = 3.5;
// For FCC the first nearest neighbor is 12
const int kNumFirstNearestNeighbors = 12;
// For FCC the first and second nearest neighbor is 18
const int kNumNearNeighbors = 18;
}// namespace FCC_const

#endif //KN_SRC_CONSTANTS_H_

