#ifndef KN_SRC_CONSTANTS_H_
#define KN_SRC_CONSTANTS_H_


// 3D simulation
const int kDimension = 3;
// three indexes
const int kXDim = 0, kYDim = 1, kZDim = 2;


namespace FCC_const {
// For FCC the first nearest neighbor is 12
const int kNumFirstNearestNeighbors = 12;
// For FCC the first and second nearest neighbor is 18
const int kNumNearNeighbors = 18;
}// namespace FCC_const

#endif //KN_SRC_CONSTANTS_H_

